#include "MHproposals_formdiss.ihs.h"

/* Shorthand. */
#define Mhead (MHp->togglehead)
#define Mtail (MHp->toggletail)

/********************
   void MH_Formation
   Propose ONLY edges not in the reference graph
***********************/
void MH_Formation (MHproposal *MHp, DegreeBound *bd, Network *nwp) 
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
    MHp->togglehead[0]=MH_FAILED;
    MHp->toggletail[0]=MH_IMPOSSIBLE;
    return;
  }

  for(trytoggle=0;trytoggle<MAX_TRIES;trytoggle++){
    /* Keep trying dyads until neither or both of nwp[0] and nwp[1] has
       the selected dyad. (That is, that dyad originally did not have an edge
       (which may have been toggled.) */
    Mtail[0] =  (nwp->directed_flag ? 1 : 2) + unif_rand() * nnodes;
    if(nwp->directed_flag){
      Mhead[0] = 1 + unif_rand() * (nnodes-1);
      if(Mhead[0]>=Mtail[0]) Mhead[0]++;
    }else{
      Mhead[0] = 1 + unif_rand() * (Mtail[0]-1);
    }
      
    if(XNOR(EdgetreeSearch(Mhead[0],Mtail[0],nwp[0].outedges),
	    EdgetreeSearch(Mhead[0],Mtail[0],nwp[1].outedges) &&
	    CheckTogglesValid(MHp, bd, nwp))) break;
  }

  /* If no valid proposal found, signal a failed proposal. */
  if(trytoggle>=MAX_TRIES) {
    Mhead[0]=MH_FAILED;
    Mtail[0]=MH_UNSUCCESSFUL;
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
void MH_FormationTNT (MHproposal *MHp, DegreeBound *bd, Network *nwp) 
{  
  Vertex head, tail;
  Edge rane, nedges, ndedges;
  double comp, odds;
  unsigned int trytoggle=0;
  static Edge ndyads;
  static Edge nnodes;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    nnodes = nwp[0].nnodes;
    ndyads = (nnodes-1)*nnodes / (nwp[0].directed_flag? 1:2);
    return;
  }
  
  nedges  = nwp[0].nedges;
  ndedges = nwp[1].nedges;

  if(nedges==ndyads && ndedges==0){  /* Attempting formation on a complete graph. */
    MHp->togglehead[0]=MH_FAILED;
    MHp->toggletail[0]=MH_IMPOSSIBLE;
    return;
  }
  
  if (ndedges > 0) {
    comp = ((double)ndedges*5.0)/((double)nedges);
    if(comp > 0.5){comp = 0.5;}
    odds = comp/(1.0-comp);
  }else{
    odds = 0.0;
  }

  for(int trytoggle = 0; trytoggle < MAX_TRIES; trytoggle++){
/*  nddyads = ndyads-nedges+ndedges; */
/*  Rprintf("comp %f nwp[0].nedges %d nwp[1].nedges %d %d %d\n",  comp, */
/*		       nwp[0].nedges, */
/*		       nwp[1].nedges, ndedges, nedges); */

    if (ndedges > 0 && unif_rand() < comp) { /* Select a new tie at random */
      rane = 1 + unif_rand() * ndedges;
      FindithEdge(&head, &tail, rane, &nwp[1]);
      /* select a dyad not in the reference network at random */
/*   Rprintf("nontie h %d t %d nwp[0].nedges %d nwp[1].nedges %d\n",  MHp->togglehead[0],   */
/*		       MHp->toggletail[0],  */
/*		       nwp[0].nedges, */
/*		       nwp[1].nedges); */
      MHp->ratio = nedges  / (odds*ndyads + nedges);
    }else{ /* select a dyad not an edge in the reference network at random */
      trytoggle=0;
      do{ /* Keep trying dyads as long as exactly one of nwp[0] and nwp[1] has
	     the selected dyad. (That is, that dyad originally had an edge (which
	     may have been toggled.) */
	head = 1 + unif_rand() * nnodes;
	while((tail = 1 + unif_rand() * nnodes) == head);
	trytoggle++;
      }while(((head > tail && !nwp->directed_flag) || 
	      XOR(EdgetreeSearch(head,tail,nwp[0].outedges),
		  EdgetreeSearch(head,tail,nwp[1].outedges))) && trytoggle<MAX_TRIES);
      
      /* Proposal failed after trying */
      if(trytoggle >= MAX_TRIES){
	MHp->togglehead[0]=MH_FAILED;
	MHp->toggletail[0]=MH_UNSUCCESSFUL; 
	return;
      }
      
      if(EdgetreeSearch(head,tail,nwp[1].outedges)!=0){
	MHp->ratio = nedges / (odds*ndyads + nedges);
      }else{
	MHp->ratio = 1.0 + (odds*ndyads)/(nedges + 1);
      }
      /*   Rprintf("nontie h %d t %d nwp[0].nedges %d nwp[1].nedges %d\n",  MHp->togglehead[0],   */
      /*		       MHp->toggletail[0],  */
      /*		       nwp[0].nedges, */
      /*		       nwp[1].nedges); */
    }

    MHp->togglehead[0]=head;
    MHp->toggletail[0]=tail;    
    if(CheckTogglesValid(MHp, bd, nwp)) break;
    
  }

  /* If tries ran out, return failure code. */
  if(trytoggle >= MAX_TRIES){
    MHp->togglehead[0]=MH_FAILED;
    MHp->toggletail[0]=MH_UNSUCCESSFUL; 
    return;
  }
}

/********************
   void MH_Dissolution
   Propose ONLY edges in the reference graph
   Any candidate edge must be either in nwp[0] or in nwp[1]. This makes
   proposing easy.
***********************/
void MH_Dissolution (MHproposal *MHp, DegreeBound *bd, Network *nwp) 
{  
  Edge rane, nedges=nwp[0].nedges, ndedges=nwp[1].nedges;
  unsigned int trytoggle;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    return;
  }
  
  if(nedges==0 && ndedges==0){ /* Attempting dissolution on an empty network. */
    Mhead[0]=MH_FAILED;
    Mtail[0]=MH_IMPOSSIBLE;
    return;
  }
  
  /* Make sure the edge is not in nwp[0] and nwp[1], which shouldn't
     happen anyway, but better safe than sorry. */
  for(trytoggle=0; trytoggle<MAX_TRIES; trytoggle++){
    rane = 1 + unif_rand() * (nedges+ndedges);
    if(rane<=nedges){
      FindithEdge(Mhead, Mtail, rane, nwp);
      if(EdgetreeSearch(Mhead[0],Mtail[0],nwp[1].outedges)==0 && CheckTogglesValid(MHp, bd, nwp)) break;
    }else{
      FindithEdge(Mhead, Mtail, rane-nedges, nwp+1);
      if(EdgetreeSearch(Mhead[0],Mtail[0],nwp[0].outedges)==0 && CheckTogglesValid(MHp, bd, nwp)) break;
    }
  }

  /* Proposal failed after trying. */
  if(trytoggle>=MAX_TRIES){ 
    Mhead[0]=MH_FAILED;
    Mtail[0]=MH_UNSUCCESSFUL;
    return;
  }

    
  MHp->ratio = 1.0;
}

/********************
   void MH_BipartiteFormation
   Propose ONLY edges not in the reference graph
***********************/
void MH_BipartiteFormation (MHproposal *MHp, DegreeBound *bd, Network *nwp) 
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
    MHp->togglehead[0]=MH_FAILED;
    MHp->toggletail[0]=MH_IMPOSSIBLE;
    return;
  }
  
  /* select a dyad not in the reference network at random */
  for(trytoggle=0;trytoggle<MAX_TRIES;trytoggle++){
    /* Keep trying dyads until either neither or both nwp[0] and nwp[1] have
	the selected dyad and the toggle satisfies degree constraints.
	(That is, that dyad originally had no edge (but	may have been toggled.) */
    Mhead[0] = 1 + unif_rand() * nb1;
    Mtail[0] = 1 + nb1 + unif_rand() * (nnodes - nb1);
    if(XNOR(EdgetreeSearch(Mhead[0],Mtail[0],nwp[0].outedges),
	    EdgetreeSearch(Mhead[0],Mtail[0],nwp[1].outedges)) &&
       CheckTogglesValid(MHp, bd, nwp)) break;
  }
  
  /* If no valid proposal found, signal a failed proposal. */
  if(trytoggle>=MAX_TRIES) {
    Mhead[0]=MH_FAILED;
    Mtail[0]=MH_UNSUCCESSFUL;
  }

  MHp->ratio = 1.0;
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
void MH_BipartiteFormationTNT (MHproposal *MHp, DegreeBound *bd, Network *nwp) 
{  
  Vertex head, tail;
  Edge rane, nedges, ndedges;
  double comp, odds;
/*  static double comp=0.99; */
/*  static double odds; */
/*  static Edge ndedges; */
  static Edge ndyads;
  static Edge nnodes;
  static Edge nb1;
  unsigned int trytoggle;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    nnodes = nwp[0].nnodes;
    nb1 = nwp[0].bipartite;
    ndyads = (nnodes-nb1)*nb1;  
    return;
  }
/*  Rprintf("nb1 %d\n",  nb1); */
  
  nedges  = nwp[0].nedges;
  ndedges = nwp[1].nedges;

  if(nwp[0].nedges==ndyads && nwp[1].nedges==0){ /* Attempting formation on a complete graph. */
    Mhead[0]=MH_FAILED;
    Mtail[0]=MH_IMPOSSIBLE;
    return;
  }

  if (ndedges > 0){
    comp = (ndedges*5.0)/(1.0*nedges);
    if(comp > 0.5){
      comp = 0.5;
    }
    odds = comp/(1.0-comp);
  }else{
    odds = 0.0;
  }
/*  nddyads = ndyads-nedges+ndedges; */
/*  Rprintf("comp %f nwp[0].nedges %d nwp[1].nedges %d %d %d\n",  comp, */
/*		       nwp[0].nedges, */
/*		       nwp[1].nedges, ndedges, nedges); */

  for(trytoggle=0; trytoggle < MAX_TRIES; trytoggle++){

    if (ndedges > 0 && unif_rand() < comp) { /* Select a new tie at random */
      rane = 1 + unif_rand() * ndedges;
      FindithEdge(MHp->togglehead, MHp->toggletail, rane, &nwp[1]);
      /* select a dyad not in the reference network at random */
/*   Rprintf("nontie h %d t %d nwp[0].nedges %d nwp[1].nedges %d\n",  MHp->togglehead[0],   */
/*		       MHp->toggletail[0],  */
/*		       nwp[0].nedges, */
/*		       nwp[1].nedges); */
    MHp->ratio = nedges  / (odds*ndyads + nedges);
    }else{ /* select a dyad not an edge in the reference network at random */
      do{ /* Keep trying dyads as long as exactly one of nwp[0] and nwp[1] has
	     the selected dyad. (That is, that dyad originally had an edge (which
	     may have been toggled.) */
	head = 1 + unif_rand() * nb1;
	tail = 1 + nb1 + unif_rand() * (nnodes - nb1);
      }while(XOR(EdgetreeSearch(head,tail,nwp[0].outedges),
		 EdgetreeSearch(head,tail,nwp[1].outedges)));
                                                                     
      MHp->togglehead[0] = head;
      MHp->toggletail[0] = tail;
      if(EdgetreeSearch(MHp->togglehead[0],MHp->toggletail[0],nwp[1].outedges)!=0){
	MHp->ratio = nedges / (odds*ndyads + nedges);
      }else{
	MHp->ratio = 1.0 + (odds*ndyads)/(nedges + 1);
      }
/*   Rprintf("nontie h %d t %d nwp[0].nedges %d nwp[1].nedges %d\n",  MHp->togglehead[0],   */
/*		       MHp->toggletail[0],  */
/*		       nwp[0].nedges, */
/*		       nwp[1].nedges); */
    }
    if(CheckTogglesValid(MHp, bd, nwp)) break;
  }

  /* If no valid proposal found, signal a failed proposal. */
  if(trytoggle>=MAX_TRIES) {
    Mhead[0]=MH_FAILED;
    Mtail[0]=MH_UNSUCCESSFUL;
  }
}

/********************
   void MH_BipartiteDissolution
   Propose ONLY edges in the reference graph
   See MH_Dissolution --- the code turns out to be identical.
***********************/
void MH_BipartiteDissolution (MHproposal *MHp, DegreeBound *bd, Network *nwp) 
{  
  MH_Dissolution(MHp,bd,nwp);
  return;
}
