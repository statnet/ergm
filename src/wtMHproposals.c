#include "wtMHproposals.h"

/* Shorthand. */
#define Mhead (MHp->togglehead)
#define Mtail (MHp->toggletail)
#define Mweight (MHp->toggleweight)

/*********************
 void MH_PseudoPoisson

 Default weighted MH algorithm for PseudoPoisson
*********************/
void MH_PseudoPoisson (WtMHproposal *MHp, WtNetwork *nwp)  {  
  Vertex head, tail;
  double oldwt, inc;
  int fvalid, trytoggle;
  
  if(MHp->ntoggles == 0) { // Initialize PseudoPoisson 
    MHp->ntoggles=1;
    return;
  }
  MHp->ratio = 1.0;
  
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

  // If we can, decrement with probability 1/2. Otherwise, increment.
  inc = (unif_rand()<0.5 || oldwt==0) ? 1 : -1;
  Mweight[0]+=inc;
  
  // Multiply the acceptance ratio by the ratio of reference measures.
  // This is y!/(y+1)! if incrementing and y!/(y-1)! if decrementing.
  MHp->ratio = inc>0 ? 1.0/Mweight[0] : oldwt;
  
}
