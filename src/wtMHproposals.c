#include "wtMHproposals.h"

/* Shorthand. */
#define Mhead (MHp->togglehead)
#define Mtail (MHp->toggletail)
#define Mweight (MHp->toggleweight)

/*********************
 void MH_Poisson

 Default MH algorithm for Poisson-reference ERGM
*********************/
void MH_Poisson(WtMHproposal *MHp, WtNetwork *nwp)  {  
  Vertex head, tail;
  double oldwt;
  int fvalid;
  
  if(MHp->ntoggles == 0) { // Initialize Poisson 
    MHp->ntoggles=1;
    return;
  }
  
  fvalid = 0;
  
  head = 1 + unif_rand() * nwp->nnodes;
  while ((tail = 1 + unif_rand() * nwp->nnodes) == head);
  if (!nwp->directed_flag && head > tail) {
    Mhead[0] = tail;
    Mtail[0] = head;
  }else{
    Mhead[0] = head;
    Mtail[0] = tail;
  }
  
  oldwt = WtGetEdge(Mhead[0],Mtail[0],nwp);

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
  Edge rane, nmissing = MHp->inputs[0];
  double oldwt;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    return;
  }

  if(nmissing==0){
    *Mhead = MH_FAILED;
    *Mtail = MH_IMPOSSIBLE;
  }


  // Note that missing edgelist is indexed from 0 but the first
  // element of MHp->inputs is the number of missing edges.
  rane = 1 + unif_rand() * nmissing;

  Mhead[0]=MHp->inputs[rane];
  Mtail[1]=MHp->inputs[nmissing+rane];

  oldwt = WtGetEdge(Mhead[0],Mtail[0],nwp);

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
  Vertex head, tail1, tail2;
  
  if(MHp->ntoggles == 0) { // Initialize Poisson 
    MHp->ntoggles=2;
    return;
  }
  
  head = 1 + unif_rand() * nwp->nnodes;
  while((tail1 = 1 + unif_rand() * nwp->nnodes) == head);
  do{
    tail2 = 1 + unif_rand() * nwp->nnodes;
  }while(tail2 == head || tail2 == tail1);
  
  Mhead[0] = Mhead[1] = head;
  Mtail[0] = tail1;
  Mtail[1] = tail2;
  
  Mweight[1] = WtGetEdge(Mhead[0],Mtail[0],nwp);
  Mweight[0] = WtGetEdge(Mhead[1],Mtail[1],nwp);
}

/*********************
 void MH_StdNormal

 Default MH algorithm for a standard-normal-reference ERGM
*********************/
void MH_StdNormal(WtMHproposal *MHp, WtNetwork *nwp)  {  
  Vertex head, tail;
  double oldwt;
  int fvalid;
  
  if(MHp->ntoggles == 0) { // Initialize Poisson 
    MHp->ntoggles=1;
    return;
  }
  
  fvalid = 0;
  
  head = 1 + unif_rand() * nwp->nnodes;
  while ((tail = 1 + unif_rand() * nwp->nnodes) == head);
  if (!nwp->directed_flag && head > tail) {
    Mhead[0] = tail;
    Mtail[0] = head;
  }else{
    Mhead[0] = head;
    Mtail[0] = tail;
  }
  
  oldwt = WtGetEdge(Mhead[0],Mtail[0],nwp);

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
  
  Mhead[0] = 1 + unif_rand() * nwp->nnodes;
  while ((Mtail[0] = 1 + unif_rand() * nwp->nnodes) == Mhead[0]);
 
  // Note that undirected networks do not make sense for this proposal.

  oldwt = WtGetEdge(Mhead[0],Mtail[0],nwp);

  // Evaluate rank constraints for this dyad.
  double ub = +HUGE_VAL, lb = -HUGE_VAL;
  int tail_class = MHp->inputs[(Mhead[0]-1)*nwp->nnodes + (Mtail[0]-1)];

  for(unsigned int t=1;t<=nwp->nnodes; t++){
    if(t==Mhead[0] || t==Mtail[0]) continue;
    unsigned int t_class=MHp->inputs[(Mhead[0]-1)*nwp->nnodes + (t-1)];

    if(t_class){
      if(t_class==tail_class+1){
	// If this alter bounds the current alter from above...
	double b = WtGetEdge(Mhead[0],t,nwp);
	if(b<ub) ub=b;
      }else if(t_class==tail_class-1){
	// If this alter bounds the current alter from below...
	double b = WtGetEdge(Mhead[0],t,nwp);
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
