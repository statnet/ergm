#include "wtMHproposals.h"
/* Shorthand. */
#define Mhead (MHp->togglehead)
#define Mtail (MHp->toggletail)
#define Mweight (MHp->toggleweight)

/*********************
 void MH_Poisson

 Default MH algorithm for Poisson-reference ERGM
*********************/
void MH_BipartitePoisson(WtMHproposal *MHp, WtNetwork *nwp)  {  
  Vertex head, tail;
  double oldwt, inc;
  int fvalid, trytoggle;
  
  if(MHp->ntoggles == 0) { // Initialize Poisson 
    MHp->ntoggles=1;
    return;
  }
  MHp->ratio = 1.0;
  
  Mhead[0] = 1 + unif_rand() * nwp->bipartite;
  Mtail[0] = 1 + nwp->bipartite + unif_rand() * (nwp->nnodes - nwp->bipartite);

  oldwt = WtGetEdge(Mhead[0],Mtail[0],nwp);

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
}

/*********************
 void MH_CompleteOrderingBipartite

 Default MH algorithm for ERGM over complete orderings
*********************/
void MH_CompleteOrderingBipartite(WtMHproposal *MHp, WtNetwork *nwp)  {  
  Vertex head, tail1, tail2;
  double oldwt, inc;
  int fvalid, trytoggle;
  
  if(MHp->ntoggles == 0) { // Initialize Poisson 
    MHp->ntoggles=2;
    return;
  }
    
  
  Mhead[0] = Mhead[1] = 1 + unif_rand() * nwp->bipartite;
  Mtail[0] = 1 + nwp->bipartite + unif_rand() * (nwp->nnodes - nwp->bipartite);
  Mtail[1] = 1 + Mtail[0] + unif_rand() * (nwp->nnodes - Mtail[0]);
  
  Mweight[1] = WtGetEdge(Mhead[0],Mtail[0],nwp);
  Mweight[0] = WtGetEdge(Mhead[1],Mtail[1],nwp);
}
