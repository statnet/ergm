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
  double oldwt, newwt;
  int fvalid, trytoggle;
  
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
  /*
  if(oldwt==0){
    // Starting with 0, can only incremement, but must adjust MH ratio.
    inc = +1;
    MHp->ratio *= 0.5;
  }else{
    // Choose whether to increment or to decrement.
    inc = (unif_rand()<0.5) ? +1 : -1;
  }

  // If going from 1 to 0, adjust the MH ratio for the forced increment.
  if(oldwt==1 && inc==-1) MHp->ratio *= 2;
  
  Mweight[0]= oldwt + inc;
  
  // Multiply the acceptance ratio by the ratio of reference measures.
  // This is y!/(y+1)! if incrementing and y!/(y-1)! if decrementing.
  MHp->ratio *= inc>0 ? 1.0/Mweight[0] : oldwt;
  */

  const double fudge = 0.5; // Mostly comes in when proposing from 0.

  do{
    newwt = rpois(oldwt + fudge);    
  }while(newwt==oldwt);
    
  //MHp->ratio *= exp((1+log(newwt + 0.5))*oldwt - (1+log(oldwt + 0.5))*newwt);
  MHp->ratio *= exp((1 + log(newwt+fudge))*oldwt - (1 + log(oldwt+fudge))*newwt) * (1-dpois(oldwt,oldwt+fudge,0))/(1-dpois(newwt,newwt+fudge,0));
  Mweight[0]=newwt;
}

/*********************
 void MH_CompleteOrdering

 Default MH algorithm for ERGM over complete orderings
*********************/
void MH_CompleteOrdering(WtMHproposal *MHp, WtNetwork *nwp)  {  
  Vertex head, tail1, tail2;
  double oldwt, inc;
  int fvalid, trytoggle;
  
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
