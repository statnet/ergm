/*
 *  File ergm/src/MHproposals_DynMoME.c
 *  Part of the statnet package, http://statnet.org
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) in
 *    http://statnet.org/attribution
 *
 *  Copyright 2012 the statnet development team
 */
#include "MHproposals_DynMoME.h"

/* Shorthand. */
#define Mtail (MHp->toggletail)
#define Mhead (MHp->togglehead)

/********************
   void MH_Formation
   Propose ONLY edges not in the reference graph
***********************/
void MH_Formation (MHproposal *MHp, Network *nwp) 
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
  
  if(nwp[0].nedges==ndyads && nwp[1].nedges==0){ /* Attempting formation on a complete graph. */
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_IMPOSSIBLE;
    return;
  }

  for(trytoggle=0;trytoggle<MAX_TRIES;trytoggle++){
    /* Keep trying dyads until neither or both of nwp[0] and nwp[1] has
       the selected dyad. (That is, that dyad originally did not have an edge
       (which may have been toggled.) */
    /* Generate. */
    Mhead[0] = 1 + unif_rand() * nnodes;
    Mtail[0] = 1 + unif_rand() * (nnodes-1);
    if(Mtail[0]>=Mhead[0]) Mtail[0]++;
    
    /* If undirected, reorder. */
    if(!nwp->directed_flag && Mhead[0]<Mtail[0]){
      Vertex tmp=Mhead[0];
      Mhead[0]=Mtail[0];
      Mtail[0]=tmp;
    }
      
    if(XNOR(EdgetreeSearch(Mtail[0],Mhead[0],nwp[0].outedges),
     EdgetreeSearch(Mtail[0],Mhead[0],nwp[1].outedges)) &&
     CheckTogglesValid(MHp, nwp)) break;
  }

  /* If no valid proposal found, signal a failed proposal. */
  if(trytoggle>=MAX_TRIES) {
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_UNSUCCESSFUL;
  }
}

/********************
   void MH_FormationTNT
   Tie/no tie:  Gives at least 50% chance of
   proposing a toggle of an existing (non-reference) edge, as opposed
   to simple random toggles that rarely do so in sparse 
   networks
   Propose ONLY edges not in the reference graph
***********************/
void MH_FormationTNT (MHproposal *MHp, Network *nwp) 
{  
  unsigned int trytoggle;
  Vertex tail, head;
  Edge nedges, ndedges, nempty;
  static double comp=0.5, odds;
  static Edge ndyads;
  static Edge nnodes;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    nnodes = nwp[0].nnodes;
    odds = comp/(1.0-comp);
    ndyads = (nnodes-1)*nnodes / (nwp[0].directed_flag? 1:2);
    return;
  }
  
  nedges  = nwp[0].nedges;
  ndedges = nwp[1].nedges;
  nempty = ndyads-nedges;

  if(nempty==0 && ndedges==0){  /* Attempting formation on a complete graph. */
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_IMPOSSIBLE;
    return;
  }
  
  for(trytoggle = 0; trytoggle < MAX_TRIES*2; trytoggle++){
    if(ndedges != 0 && unif_rand() < comp) { /* Select a discordant dyad at random */
      GetRandEdge(&tail, &head, nwp+1);
      
      MHp->logratio += log(ndedges / ((double)nempty + 1)/comp * ((ndedges==1)? 1 : (1-comp)));
    }else{ /* select an empty dyad in nwp[0] at random */
      do{ /* Keep trying dyads as long as it's an edge in nwp[0]. */
        head = 1 + unif_rand() * nnodes;
        tail = 1 + unif_rand() * (nnodes-1);
        if(tail>=head) tail++;
        trytoggle++;
      }while(((tail > head && !nwp->directed_flag) || 
          EdgetreeSearch(tail,head,nwp[0].outedges)) && trytoggle<MAX_TRIES);
      
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
//   Rprintf("nedges %d reference nddyads %d h %d t %d MHp->ratio %f\n", 
//	    nwp[0].nedges, nwp[1].nedges, tail, head, MHp->ratio); 
  
  /* If tries ran out, return failure code. */
  if(trytoggle >= MAX_TRIES*2){
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_UNSUCCESSFUL; 
    return;
  }
}

/********************
   void MH_Dissolution
   Propose ONLY edges in the reference graph
   Any candidate edge must be either in nwp[0] or in nwp[1]. This makes
   proposing easy.
***********************/
void MH_Dissolution (MHproposal *MHp, Network *nwp) 
{  
  Edge nedges=nwp[0].nedges, ndedges=nwp[1].nedges;
  unsigned int trytoggle;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    return;
  }
  
  if(nedges==0 && ndedges==0){ /* Attempting dissolution on an empty network. */
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_IMPOSSIBLE;
    return;
  }
  
  /* Make sure the edge is not in nwp[0] and nwp[1], which shouldn't
     happen anyway, but better safe than sorry. */
  for(trytoggle=0; trytoggle<MAX_TRIES; trytoggle++){
    Edge rane = 1 + unif_rand() * (nedges+ndedges);
    if(rane<=nedges){
      GetRandEdge(Mtail, Mhead, nwp);
      if(EdgetreeSearch(Mtail[0],Mhead[0],nwp[1].outedges)==0 && CheckTogglesValid(MHp, nwp)) break;
    }else{
      GetRandEdge(Mtail, Mhead, nwp+1);
      if(EdgetreeSearch(Mtail[0],Mhead[0],nwp[0].outedges)==0 && CheckTogglesValid(MHp, nwp)) break;
    }
  }

  /* Proposal failed after trying. */
  if(trytoggle>=MAX_TRIES){ 
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_UNSUCCESSFUL;
    return;
  }
}

/********************
   void MH_BipartiteFormation
   Propose ONLY edges not in the reference graph
***********************/
void MH_BipartiteFormation (MHproposal *MHp, Network *nwp) 
{  
  static Edge nnodes;
  static Edge nb1;
  unsigned int trytoggle;
  static Edge ndyads;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    nnodes = nwp[0].nnodes;
    nb1 = nwp[0].bipartite;
    ndyads = (nnodes-nb1)*nb1;
    return;
  }

  if(nwp[0].nedges==ndyads && nwp[1].nedges==0){ /* Attempting formation on a complete graph. */
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_IMPOSSIBLE;
    return;
  }
  
  /* select a dyad not in the reference network at random */
  for(trytoggle=0;trytoggle<MAX_TRIES;trytoggle++){
    /* Keep trying dyads until either neither or both nwp[0] and nwp[1] have
	the selected dyad and the toggle satisfies degree constraints.
	(That is, that dyad originally had no edge (but	may have been toggled.) */
    Mtail[0] = 1 + unif_rand() * nb1;
    Mhead[0] = 1 + nb1 + unif_rand() * (nnodes - nb1);
    if(XNOR(EdgetreeSearch(Mtail[0],Mhead[0],nwp[0].outedges),
	    EdgetreeSearch(Mtail[0],Mhead[0],nwp[1].outedges)) &&
       CheckTogglesValid(MHp, nwp)) break;
  }
  
  /* If no valid proposal found, signal a failed proposal. */
  if(trytoggle>=MAX_TRIES) {
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_UNSUCCESSFUL;
  }

/*   Rprintf("reference nddyads %d MHp->ratio %f\n", */
/*	    nwp[1].nedges, MHp->ratio); */
}

/********************
   void MH_BipartiteFormationTNT
   Tie/no tie:  Gives at least 50% chance of
   proposing a toggle of an existing (non-reference) edge, as opposed
   to simple random toggles that rarely do so in sparse 
   networks
   Propose ONLY edges not in the reference graph
***********************/
void MH_BipartiteFormationTNT (MHproposal *MHp, Network *nwp) 
{  
  unsigned int trytoggle;
  Vertex tail, head;
  Edge nedges, ndedges, nempty;
  static double comp=0.5, odds;
  static Edge ndyads;
  static Edge nnodes;
  static Edge nb1;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    nnodes = nwp[0].nnodes;
    odds = comp/(1.0-comp);
    nb1 = nwp[0].bipartite;
    ndyads = (nnodes-nb1)*nb1;  
    return;
  }
  
  nedges  = nwp[0].nedges;
  ndedges = nwp[1].nedges;
  nempty = ndyads-nedges;
  
  if(nempty==0 && ndedges==0){ /* Attempting formation on a complete graph. */
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_IMPOSSIBLE;
    return;
  }
  
  for(trytoggle=0; trytoggle < MAX_TRIES*2; trytoggle++){
    if(ndedges > 0 && unif_rand() < comp) { /* Select a discordant dyad at random */
      GetRandEdge(&tail, &head, nwp+1);
      
      MHp->logratio += log(ndedges  / ((double)nempty+1)/comp *((ndedges==1)? 1 : (1-comp)));
    }else{ /* select an empty dyad in nwp[0] at random */
      do{ /* Keep trying dyads as long as it's an edge in nwp[0]. */
        tail = 1 + unif_rand() * nb1;
        head = 1 + nb1 + unif_rand() * (nnodes - nb1);
        trytoggle++;
      }while(EdgetreeSearch(tail,head,nwp[0].outedges) && trytoggle<MAX_TRIES);
      
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
    
    Mtail[0] = tail;
    Mhead[0] = head;
    
    if(CheckTogglesValid(MHp, nwp)) break;
  }
  
  /* If no valid proposal found, signal a failed proposal. */
  if(trytoggle>=MAX_TRIES*2) {
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_UNSUCCESSFUL;
  }
}

/********************
   void MH_BipartiteDissolution
   Propose ONLY edges in the reference graph
   See MH_Dissolution --- the code turns out to be identical.
***********************/
void MH_BipartiteDissolution (MHproposal *MHp, Network *nwp) 
{  
  MH_Dissolution(MHp, nwp);
  return;
}
