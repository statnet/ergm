#include "wtMHproposals.h"
/* Shorthand. */
#define Mhead (MHp->togglehead)
#define Mtail (MHp->toggletail)
#define Mweight (MHp->toggleweight)

/*********************
 void MH_Poisson

 Default MH algorithm for Poisson-reference ERGM
*********************/
void MH_BipartitePoisson (WtMHproposal *MHp, WtNetwork *nwp)  {  
  Vertex head, tail;
  double oldwt, inc;
  int fvalid, trytoggle;
  
  if(MHp->ntoggles == 0) { // Initialize Poisson 
    MHp->ntoggles=1;
    return;
  }
  MHp->ratio = 1.0;
  
  head = 1 + unif_rand() * nwp->nnodes;
  Mhead[0] = 1 + unif_rand() * nwp->bipartite;
  Mtail[0] = 1 + nwp->bipartite + 
    unif_rand() * (nwp->nnodes - nwp->bipartite);

  oldwt = WtGetEdge(Mhead[0],Mtail[0],nwp);

  // If we can, decrement with probability 1/2. Otherwise, increment.
  inc = (unif_rand()<0.5) ? +1 : -1;

  Mweight[0]=oldwt + inc;
  
  // If the proposed weight is outside of the support, automatically
  // reject the proposal.
  if(Mweight[0]<0){
    MHp->ratio=0;
  }
  else{
    // Multiply the acceptance ratio by the ratio of reference measures.
  // This is y!/(y+1)! if incrementing and y!/(y-1)! if decrementing.
    MHp->ratio = inc>0 ? 1.0/Mweight[0] : oldwt;
  }
  
}
