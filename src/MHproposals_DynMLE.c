/*
 *  File ergm/src/MHproposals_DynMLE.c
 *  Part of the statnet package, http://statnet.org
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) in
 *    http://statnet.org/attribution
 *
 *  Copyright 2012 the statnet development team
 */
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
  static Edge nedges0;

  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    nnodes = nwp[0].nnodes;
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

/********************
   void MH_FormationMLE
   Propose ONLY edges not in the reference graph
***********************/
void MH_FormationMLETNT(MHproposal *MHp, Network *nwp) 
{  
  static Vertex nnodes;
  unsigned int trytoggle;
  Vertex tail,head;
  static Edge ndyads;
  static double comp=0.5, odds;
  static Network discord;

  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    nnodes = nwp->nnodes;
    odds = comp/(1.0-comp);
    ndyads = (nnodes-1)*nnodes / (nwp->directed_flag? 1:2);

    Edge nedges0 = MHp->inputs[0];
    MHp->discord = (Network**) calloc(2,sizeof(Network*)); // A space for the sentinel NULL pointer.
    MHp->discord[0] = &discord;
    discord = NetworkInitializeD(MHp->inputs+1, MHp->inputs+1+nedges0, nedges0, nnodes, nwp->directed_flag, nwp->bipartite, 0, 0, NULL);
   
    for(Edge i=0; i<nwp->nedges; i++){
      FindithEdge(&tail, &head, i+1, nwp);
      ToggleEdge(tail,head,&discord);
    }
    
    return;
  }
  

  Edge nedges = nwp->nedges, ndedges = discord.nedges, nempty = ndyads-nedges;

  if(nempty==0 && ndedges==0){ /* Attempting formation on a complete graph. */
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_IMPOSSIBLE;
    return;
  }

  for(trytoggle=0;trytoggle<MAX_TRIES*2;trytoggle++){
    if(ndedges != 0 && unif_rand() < comp) { /* Select a discordant dyad at random */
      GetRandEdge(&tail, &head, &discord);
      
      MHp->logratio += log(ndedges / ((double)nempty + 1)/comp * ((ndedges==1)? 1 : (1-comp)));
    }else{ /* select an empty dyad in nwp[0] at random */
      do{ /* Keep trying dyads as long as it's an edge in nwp */
        head = 1 + unif_rand() * nnodes;
        tail = 1 + unif_rand() * (nnodes-1);
        if(tail>=head) tail++;
        trytoggle++;
      }while(((tail > head && !nwp->directed_flag) || 
	      EdgetreeSearch(tail,head,nwp->outedges)) && trytoggle<MAX_TRIES);
      
      /* Proposal failed after trying */
      if(trytoggle >= MAX_TRIES*2){
        Mtail[0]=MH_FAILED;
        Mhead[0]=MH_UNSUCCESSFUL; 
        return;
      }
      
      if(ndedges==0){
        MHp->logratio += log(nempty*comp);
      }else{
        MHp->logratio += log(((double)nempty)/(ndedges+1) *odds);
      }
    }
    
    Mtail[0]=tail;
    Mhead[0]=head;
    
    if(CheckTogglesValid(MHp, nwp)) break;
  }

  /* If no valid proposal found, signal a failed proposal. */
  if(trytoggle>=MAX_TRIES*2) {
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_UNSUCCESSFUL;
  }
}
