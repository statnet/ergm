#include "MHproposals_DynMLE.h"
#include "edgelist.h"

/* Shorthand. */
#define Mtail (MHp->toggletail)
#define Mhead (MHp->togglehead)

/********************
   void MH_FormationMLE
   Propose ONLY edges not in the reference graph
***********************/
void MH_FormationMLE (MHproposal *MHp, Network *nwp) 
{  
  static Vertex nnodes;
  unsigned int trytoggle;
  static Edge ndyads, nedges0;

  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    nnodes = nwp[0].nnodes;
    ndyads = (nnodes-1)*nnodes / (nwp[0].directed_flag? 1:2);
    nedges0 = MHp->inputs[0];
    return;
  }
  
  if(nedges0==ndyads){ /* Attempting formation on a complete graph. */
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_IMPOSSIBLE;
    return;
  }

  for(trytoggle=0;trytoggle<MAX_TRIES;trytoggle++){
    /* Keep trying dyads until a one that is not an edge in the reference network is found. */
    /* Generate. */
    do{
      Mhead[0] = 1 + unif_rand() * nnodes;
      Mtail[0] = 1 + unif_rand() * nnodes;
    }while(Mtail[0]==Mhead[0]);
    
    /* If undirected, reorder. */
    if(!nwp->directed_flag && Mhead[0]<Mtail[0]){
      Vertex tmp=Mhead[0];
      Mhead[0]=Mtail[0];
      Mtail[0]=tmp;
    }
      
    if(!dEdgeListSearch(Mtail[0],Mhead[0],MHp->inputs) &&
	    CheckTogglesValid(MHp, nwp)) break;
  }

  /* If no valid proposal found, signal a failed proposal. */
  if(trytoggle>=MAX_TRIES) {
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_UNSUCCESSFUL;
  }
}

/********************
   void MH_DissolutionMLE
   Propose ONLY edges not in the reference graph
***********************/
void MH_DissolutionMLE (MHproposal *MHp, Network *nwp) 
{  
  static Vertex nnodes;
  unsigned int trytoggle;
  static Edge ndyads, nedges0;

  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    nnodes = nwp[0].nnodes;
    ndyads = (nnodes-1)*nnodes / (nwp[0].directed_flag? 1:2);
    nedges0 = MHp->inputs[0];
    return;
  }
  
  if(nedges0==0){ /* Attempting dissolution on a complete graph. */
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_IMPOSSIBLE;
    return;
  }

  for(trytoggle=0;trytoggle<MAX_TRIES;trytoggle++){
    /* Select a dyad at random that is in the reference graph. (We
       have a convenient sampling frame.) */
    /* Generate. */
    Edge rane = 1 + unif_rand() * nedges0;
    Mtail[0]=MHp->inputs[rane];
    Mhead[0]=MHp->inputs[nedges0+rane];
    
    if(CheckTogglesValid(MHp, nwp)) break;
  }

  /* If no valid proposal found, signal a failed proposal. */
  if(trytoggle>=MAX_TRIES) {
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_UNSUCCESSFUL;
  }
}
