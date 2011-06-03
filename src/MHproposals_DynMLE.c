#include "MHproposals_DynMLE.h"

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
  static Edge ndyads;

  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    nnodes = nwp[0].nnodes;
    ndyads = (nnodes-1)*nnodes / (nwp[0].directed_flag? 1:2);
    return;
  }
  
  if(nwp[2].nedges==ndyads){ /* Attempting formation on a complete graph. */
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_IMPOSSIBLE;
    return;
  }

  for(trytoggle=0;trytoggle<MAX_TRIES;trytoggle++){
    /* Keep trying dyads until neither or both of nwp[0] and nwp[2] has
       the selected dyad. (That is, that dyad originally did not have an edge
       (which may have been toggled.) */
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
      
    if((!EdgetreeSearch(Mtail[0],Mhead[0],nwp[0].outedges)|
        !EdgetreeSearch(Mtail[0],Mhead[0],nwp[2].outedges)) &&
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
  static Edge ndyads;

  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    nnodes = nwp[0].nnodes;
    ndyads = (nnodes-1)*nnodes / (nwp[0].directed_flag? 1:2);
    return;
  }
  
  if(nwp[2].nedges==0){ /* Attempting dissolution on a complete graph. */
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_IMPOSSIBLE;
    return;
  }

  for(trytoggle=0;trytoggle<MAX_TRIES;trytoggle++){
    /* Keep trying dyads until neither or both of nwp[0] and nwp[2] has
       the selected dyad. (That is, that dyad originally did not have an edge
       (which may have been toggled.) */
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
      
    if(( EdgetreeSearch(Mtail[0],Mhead[0],nwp[0].outedges)|
         EdgetreeSearch(Mtail[0],Mhead[0],nwp[2].outedges)) &&
	    CheckTogglesValid(MHp, nwp)) break;
  }

  /* If no valid proposal found, signal a failed proposal. */
  if(trytoggle>=MAX_TRIES) {
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_UNSUCCESSFUL;
  }
}

/********************
   void MH_FormationNonObservedMLE
   Propose ONLY edges not in the reference graph
***********************/
void MH_FormationNonObservedMLE (MHproposal *MHp, Network *nwp) 
{  
  static Vertex nnodes;
  unsigned int trytoggle;
  static Edge ndyads;
  Edge nmissing = MHp->inputs[0];

  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    nnodes = nwp[0].nnodes;
    ndyads = (nnodes-1)*nnodes / (nwp[0].directed_flag? 1:2);
    if(nmissing==0){
      *MHp->toggletail = MH_FAILED;
      *MHp->togglehead = MH_IMPOSSIBLE;
    }
    return;
  }
  
  if(nwp[2].nedges==ndyads){ /* Attempting formation on a complete graph. */
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_IMPOSSIBLE;
    return;
  }

  for(trytoggle=0;trytoggle<MAX_TRIES;trytoggle++){
    /* Keep trying dyads until neither or both of nwp[0] and nwp[2] has
       the selected dyad. (That is, that dyad originally did not have an edge
       (which may have been toggled.) */
    /* Generate. */
  
    // Note that missing edgelist is indexed from 0 but the first
    // element of MHp->inputs is the number of missing edges.
    Edge rane = 1 + unif_rand() * nmissing;
    Mtail[0]=MHp->inputs[rane];
    Mhead[0]=MHp->inputs[nmissing+rane];
      
    if((!EdgetreeSearch(Mtail[0],Mhead[0],nwp[0].outedges)||
        !EdgetreeSearch(Mtail[0],Mhead[0],nwp[2].outedges)) &&
	    CheckTogglesValid(MHp, nwp)) break;
  }

  /* If no valid proposal found, signal a failed proposal. */
  if(trytoggle>=MAX_TRIES) {
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_UNSUCCESSFUL;
  }
}

/********************
   void MH_DissolutionNonObservedMLE
   Propose ONLY edges not in the reference graph
***********************/
void MH_DissolutionNonObservedMLE (MHproposal *MHp, Network *nwp) 
{  
  static Vertex nnodes;
  unsigned int trytoggle;
  static Edge ndyads;
  Edge nmissing=MHp->inputs[0];

  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    nnodes = nwp[0].nnodes;
    ndyads = (nnodes-1)*nnodes / (nwp[0].directed_flag? 1:2);
    if(nmissing==0){
      *MHp->toggletail = MH_FAILED;
      *MHp->togglehead = MH_IMPOSSIBLE;
    }
    return;
  }
  
  if(nwp[2].nedges==0){ /* Attempting dissolution on a empty graph. */
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_IMPOSSIBLE;
    return;
  }

  for(trytoggle=0;trytoggle<MAX_TRIES;trytoggle++){
    /* Keep trying dyads until neither or both of nwp[0] and nwp[2] has
       the selected dyad. (That is, that dyad originally did not have an edge
       (which may have been toggled.) */
    /* Generate. */
  
    // Note that missing edgelist is indexed from 0 but the first
    // element of MHp->inputs is the number of missing edges.
    Edge rane = 1 + unif_rand() * nmissing;
    Mtail[0]=MHp->inputs[rane];
    Mhead[0]=MHp->inputs[nmissing+rane];
      
    if(( EdgetreeSearch(Mtail[0],Mhead[0],nwp[0].outedges)||
         EdgetreeSearch(Mtail[0],Mhead[0],nwp[2].outedges)) &&
	    CheckTogglesValid(MHp, nwp)) break;
  }

  /* If no valid proposal found, signal a failed proposal. */
  if(trytoggle>=MAX_TRIES) {
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_UNSUCCESSFUL;
  }
}
