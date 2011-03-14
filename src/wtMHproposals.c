#include "wtMHproposals.h"

/* Shorthand. */
#define Mtail (MHp->toggletail)
#define Mhead (MHp->togglehead)
#define Mweight (MHp->toggleweight)

/*********************
 void MH_Poisson

 Default MH algorithm for Poisson-reference ERGM
*********************/
void MH_Poisson(WtMHproposal *MHp, WtNetwork *nwp)  {  
  Vertex tail, head;
  double oldwt;
  int fvalid;
  
  if(MHp->ntoggles == 0) { // Initialize Poisson 
    MHp->ntoggles=1;
    return;
  }
  
  fvalid = 0;
  
  tail = 1 + unif_rand() * nwp->nnodes;
  while ((head = 1 + unif_rand() * nwp->nnodes) == tail);
  if (!nwp->directed_flag && tail > head) {
    Mtail[0] = head;
    Mhead[0] = tail;
  }else{
    Mtail[0] = tail;
    Mhead[0] = head;
  }
  
  oldwt = WtGetEdge(Mtail[0],Mhead[0],nwp);

  const double fudge = 0.5; // Mostly comes in when proposing from 0.

  do{
    Mweight[0] = rpois(oldwt + fudge);    
  }while(Mweight[0]==oldwt);
    
  MHp->logratio += (1 + log(Mweight[0]+fudge))*oldwt - (1 + log(oldwt+fudge))*Mweight[0] + log(1-dpois(oldwt,oldwt+fudge,0)) - log(1-dpois(Mweight[0],Mweight[0]+fudge,0));
}

/*********************
 void MH_PoissonNonObserved

 Missing data MH algorithm for Poisson-reference ERGM on bipartite networks.
*********************/
void MH_PoissonNonObserved(WtMHproposal *MHp, WtNetwork *nwp)  {  
  Edge nmissing = MHp->inputs[0];

  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    return;
  }

  if(nmissing==0){
    *Mtail = MH_FAILED;
    *Mhead = MH_IMPOSSIBLE;
  }


  // Note that missing edgelist is indexed from 0 but the first
  // element of MHp->inputs is the number of missing edges.
  Edge rane = 1 + unif_rand() * nmissing;

  Mtail[0]=MHp->inputs[rane];
  Mhead[1]=MHp->inputs[nmissing+rane];

  double oldwt = WtGetEdge(Mtail[0],Mhead[0],nwp);

  const double fudge = 0.5; // Mostly comes in when proposing from 0.

  do{
    Mweight[0] = rpois(oldwt + fudge);    
  }while(Mweight[0]==oldwt);
    
  MHp->logratio += (1 + log(Mweight[0]+fudge))*oldwt - (1 + log(oldwt+fudge))*Mweight[0] + log(1-dpois(oldwt,oldwt+fudge,0)) - log(1-dpois(Mweight[0],Mweight[0]+fudge,0));
}


/*********************
 void MH_CompleteOrdering

 Default MH algorithm for ERGM over complete orderings
*********************/
void MH_CompleteOrdering(WtMHproposal *MHp, WtNetwork *nwp)  {  
  Vertex tail, head1, head2;
  
  if(MHp->ntoggles == 0) { // Initialize Poisson 
    MHp->ntoggles=2;
    return;
  }
  
  tail = 1 + unif_rand() * nwp->nnodes;
  while((head1 = 1 + unif_rand() * nwp->nnodes) == tail);
  do{
    head2 = 1 + unif_rand() * nwp->nnodes;
  }while(head2 == tail || head2 == head1);
  
  Mtail[0] = Mtail[1] = tail;
  Mhead[0] = head1;
  Mhead[1] = head2;
  
  Mweight[1] = WtGetEdge(Mtail[0],Mhead[0],nwp);
  Mweight[0] = WtGetEdge(Mtail[1],Mhead[1],nwp);
}

/*********************
 void MH_StdNormal

 Default MH algorithm for a standard-normal-reference ERGM
*********************/
void MH_StdNormal(WtMHproposal *MHp, WtNetwork *nwp)  {  
  Vertex tail, head;
  double oldwt;
  int fvalid;
  
  if(MHp->ntoggles == 0) { // Initialize Poisson 
    MHp->ntoggles=1;
    return;
  }
  
  fvalid = 0;
  
  tail = 1 + unif_rand() * nwp->nnodes;
  while ((head = 1 + unif_rand() * nwp->nnodes) == tail);
  if (!nwp->directed_flag && tail > head) {
    Mtail[0] = head;
    Mhead[0] = tail;
  }else{
    Mtail[0] = tail;
    Mhead[0] = head;
  }
  
  oldwt = WtGetEdge(Mtail[0],Mhead[0],nwp);

  const double propsd = 0.2; // This ought to be tunable.

  Mweight[0] = rnorm(oldwt, propsd);    
  
  // Symmetric proposal, but depends on the reference measure
  MHp->logratio += -(Mweight[0]*Mweight[0]-oldwt*oldwt)/2;
}

/*********************
 void MH_StdNormalRank

 Default MH algorithm for a standard-normal-reference ERGM
 subject to the constraint that the ranking of alters for
 each ego is preserved.
*********************/
void MH_StdNormalRank(WtMHproposal *MHp, WtNetwork *nwp)  {  
  double oldwt;
  int fvalid;
  
  if(MHp->ntoggles == 0) { // Initialize StdNormalRank
    MHp->ntoggles=1;
    return;
  }
  
  fvalid = 0;
  
  Mtail[0] = 1 + unif_rand() * nwp->nnodes;
  while ((Mhead[0] = 1 + unif_rand() * nwp->nnodes) == Mtail[0]);
 
  // Note that undirected networks do not make sense for this proposal.

  oldwt = WtGetEdge(Mtail[0],Mhead[0],nwp);

  // Evaluate rank constraints for this dyad.
  double ub = +HUGE_VAL, lb = -HUGE_VAL;
  int head_class = MHp->inputs[(Mtail[0]-1)*nwp->nnodes + (Mhead[0]-1)];

  for(unsigned int h=1;h<=nwp->nnodes; h++){
    if(h==Mtail[0] || h==Mhead[0]) continue;
    unsigned int h_class=MHp->inputs[(Mtail[0]-1)*nwp->nnodes + (h-1)];

    if(h_class){
      if(h_class==head_class+1){
	// If this alter bounds the current alter from above...
	double b = WtGetEdge(Mtail[0],h,nwp);
	if(b<ub) ub=b;
      }else if(h_class==head_class-1){
	// If this alter bounds the current alter from below...
	double b = WtGetEdge(Mtail[0],h,nwp);
	if(b>lb) lb=b;
      }
    }
  }

  const double propsd = 1; // This ought to be tunable.
  // Note that for a given dyad, which branch gets taken here is always the same.
  if(ub==+HUGE_VAL && lb==-HUGE_VAL){
    // No constraint -> a normal jump
    Mweight[0] = rnorm(oldwt, propsd);
  }else if(ub==+HUGE_VAL){
    // Only lower bound -> a constrained normal jump
    do{ Mweight[0] = rnorm(oldwt, propsd); }while(Mweight[0]<lb);
      
    MHp->logratio += dnorm(oldwt, Mweight[0], propsd, 1) - pnorm(lb, Mweight[0], propsd, 0, 1);
    MHp->logratio -= dnorm(Mweight[0], oldwt, propsd, 1) - pnorm(lb, oldwt, propsd, 0, 1);
  }else if(lb==-HUGE_VAL){
    // Only upper bound -> a constrained normal jump
    do{ Mweight[0] = rnorm(oldwt, propsd); }while(Mweight[0]>ub);
      
    MHp->logratio += dnorm(oldwt, Mweight[0], propsd, 1) - pnorm(ub, Mweight[0], propsd, 1, 1);
    MHp->logratio -= dnorm(Mweight[0], oldwt, propsd, 1) - pnorm(ub, oldwt, propsd, 1, 1);
  }else{
    // Bounded from both sides -> uniform
    Mweight[0] = runif(lb,ub);
  }

  // Depends on the reference measure
  MHp->logratio += -(Mweight[0]*Mweight[0]-oldwt*oldwt)/2;
}
