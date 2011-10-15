#include "wtMHproposals.h"
/* Shorthand. */
#define Mtail (MHp->toggletail)
#define Mhead (MHp->togglehead)
#define Mweight (MHp->toggleweight)

/*********************
 void MH_BipartitePoisson

 Default MH algorithm for Poisson-reference ERGM
*********************/
void MH_BipartitePoisson(WtMHproposal *MHp, WtNetwork *nwp)  {  
  double oldwt;
  
  if(MHp->ntoggles == 0) { // Initialize Poisson 
    MHp->ntoggles=1;
    return;
  }
  
  Mtail[0] = 1 + unif_rand() * nwp->bipartite;
  Mhead[0] = 1 + nwp->bipartite + unif_rand() * (nwp->nnodes - nwp->bipartite);

  oldwt = WtGetEdge(Mtail[0],Mhead[0],nwp);

  const double fudge = 0.5; // Mostly comes in when proposing from 0.

  do{
    Mweight[0] = rpois(oldwt + fudge);    
  }while(Mweight[0]==oldwt);
    
  MHp->logratio += (1 + log(Mweight[0]+fudge))*oldwt - (1 + log(oldwt+fudge))*Mweight[0] + log(1-dpois(oldwt,oldwt+fudge,0)) - log(1-dpois(Mweight[0],Mweight[0]+fudge,0));
}

/*********************
 void MH_BipartitePoissonNonObserved

 Missing data MH algorithm for Poisson-reference ERGM on bipartite networks.
 Completely identical to the non-bipartite version.
*********************/
void MH_BipartitePoissonNonObserved(WtMHproposal *MHp, WtNetwork *nwp){ MH_PoissonNonObserved(MHp, nwp); }

/*********************
 void MH_CompleteOrderingBipartite

 Default MH algorithm for ERGM over complete orderings
*********************/
void MH_CompleteOrderingBipartite(WtMHproposal *MHp, WtNetwork *nwp)  {  
  
  if(MHp->ntoggles == 0) { // Initialize Poisson 
    MHp->ntoggles=2;
    return;
  }
    
  
  Mtail[0] = Mtail[1] = 1 + unif_rand() * nwp->bipartite;
  Mhead[0] = 1 + nwp->bipartite + unif_rand() * (nwp->nnodes - nwp->bipartite);
  Mhead[1] = 1 + Mhead[0] + unif_rand() * (nwp->nnodes - Mhead[0]);
  
  Mweight[1] = WtGetEdge(Mtail[0],Mhead[0],nwp);
  Mweight[0] = WtGetEdge(Mtail[1],Mhead[1],nwp);
}

/*********************
 void MH_CompleteOrderingEquivalentBipartite

Completely identical to non-bipartite
*********************/

void MH_CompleteOrderingEquivalentBipartite(WtMHproposal *MHp, WtNetwork *nwp)  {  
  MH_CompleteOrderingEquivalent(MHp, nwp);
}

/*********************
 void MH_BipartiteStdNormal

 Default MH algorithm for standard-normal-reference ERGM
*********************/
void MH_BipartiteStdNormal(WtMHproposal *MHp, WtNetwork *nwp)  {  
  double oldwt;
  
  if(MHp->ntoggles == 0) { // Initialize Poisson 
    MHp->ntoggles=1;
    return;
  }
  
  Mtail[0] = 1 + unif_rand() * nwp->bipartite;
  Mhead[0] = 1 + nwp->bipartite + unif_rand() * (nwp->nnodes - nwp->bipartite);

  oldwt = WtGetEdge(Mtail[0],Mhead[0],nwp);

  const double propsd = 0.2;

  Mweight[0] = rnorm(oldwt, propsd); // This ought to be tunable.  
  
  // Symmetric proposal, but depends on the reference measure
  MHp->logratio += -(Mweight[0]*Mweight[0]-oldwt*oldwt)/2;
}
