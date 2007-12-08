#include "MHproposals_bipartite.h" 
#include "MHproposals.h" 

/* Shorthand. */
#define Mhead (MHp->togglehead)
#define Mtail (MHp->toggletail)


/*********************
 void MH_bipartite

 Default MH algorithm for bipartite networks
*********************/
void MH_Bipartiterandomtoggle (MHproposal *MHp, DegreeBound *bd, Network *nwp)  {  
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    return;
  }
  MHp->ratio = 1.0;
  MHp->togglehead[0] = 1 + unif_rand() * nwp->bipartite;
  MHp->toggletail[0] = 1 + nwp->bipartite + 
                       unif_rand() * (nwp->nnodes - nwp->bipartite);
}

/********************
   void MH_BipartiteConstantEdges
   Chooses a pair of toggles - one a tie and one not. 
***********************/
void MH_BipartiteConstantEdges (MHproposal *MHp, DegreeBound *bd, Network *nwp)  {  
  Vertex head, tail;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=2;    
    return;
  } /* Note:  This proposal cannot be used for full or empty observed graphs.
       If desired, we could check for this at initialization phase. 
       (For now, however, no way to easily return an error message and stop.)*/
  MHp->ratio=1.0;   
  /* First, select edge at random */
  FindithEdge(MHp->togglehead, MHp->toggletail, 1+(nwp->nedges)*unif_rand(), nwp);
  /* Second, select dyad at random until it has no edge */

  do{
/*    head = 1 + unif_rand() * nwp->nnodes; */
/*    tail = 1 + unif_rand() * nwp->nnodes; */
    head = 1 + unif_rand() * nwp->bipartite;
    tail = 1 + nwp->bipartite + 
      unif_rand() * (nwp->nnodes - nwp->bipartite);
  }while(EdgetreeSearch(head, tail, nwp->outedges) != 0);
      
  MHp->togglehead[1]=head;
  MHp->toggletail[1]=tail;
}

/********************
   void MH_BipartiteTNT
   Tie/no tie:  Gives at least 50% chance of
   proposing a toggle of an existing edge, as opposed
   to simple random toggles that rarely do so in sparse 
   networks
***********************/
void MH_BipartiteTNT (MHproposal *MHp, DegreeBound *bd, Network *nwp) 
{  
  Vertex head, tail;
  Edge rane, nedges=nwp->nedges;
  unsigned trytoggle;
  static double comp=0.5;
  static double odds;
  static Edge ndyads;
  static Edge nnodes;
  static Edge nb1;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    odds = comp/(1.0-comp);
    nnodes = nwp[0].nnodes;
    nb1 = nwp[0].bipartite;
    ndyads = (nnodes-nb1)*nb1;  
    return;
  }
  
  for(trytoggle = 0; trytoggle < MAX_TRIES; trytoggle++){
    if (unif_rand() < comp && nedges > 0) { /* Select a tie at random */
      rane = 1 + unif_rand() * nedges;
      FindithEdge(MHp->togglehead, MHp->toggletail, rane, nwp);
      MHp->ratio = nedges  / (odds*ndyads + nedges);
    }else{ /* Select a dyad at random */
      head = 1 + unif_rand() * nb1;
      tail = 1 + nb1 + unif_rand() * (nnodes - nb1);
      MHp->togglehead[0] = head;
      MHp->toggletail[0] = tail;
      if(EdgetreeSearch(MHp->togglehead[0],MHp->toggletail[0],nwp->outedges)!=0){
	MHp->ratio = nedges / (odds*ndyads + nedges);
      }else{
	MHp->ratio = 1.0 + (odds*ndyads)/(nedges + 1);
      }
    }
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
   void MH_BipartiteHammingConstantEdges
   Chooses a pair of toggles - one a tie and one not. 
   MSH: The name Hamming is a hack for the Hamming proposals
        It is no different the MH_BipartiteConstantEdges
***********************/
void MH_BipartiteHammingConstantEdges (MHproposal *MHp, DegreeBound *bd, Network *nwp) 
{  
  Vertex head, tail;
  Edge rane, nedges=nwp[0].nedges, nddyads=nwp[1].nedges;
  int nde, ndn, nce, ncn;
  static double comp=0.5;
  static double odds;
  static Edge ndyads;
  static Edge nnodes;
  static Edge nb1;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=2;
    odds = comp/(1.0-comp);
    nnodes = nwp[0].nnodes;
    nb1 = nwp[0].bipartite;
    ndyads = (nnodes-nb1)*nb1;  
    return;
  }
  
  if (unif_rand() < comp && nddyads > 0) { /* Select a discordant pair of tie/nontie at random */
    /* First, select discord edge at random */
    do{
      rane = 1 + unif_rand() * nddyads;
      FindithEdge(MHp->togglehead, MHp->toggletail, rane, &nwp[1]);
    }while(EdgetreeSearch(MHp->togglehead[0], MHp->toggletail[0], nwp[0].outedges) == 0);

    head=MHp->togglehead[0];
    tail=MHp->toggletail[0];
    /* Next, select discord non-edge at random */

    do{
      rane = 1 + unif_rand() * nddyads;
      FindithEdge(MHp->togglehead, MHp->toggletail, rane, &nwp[1]);
    }while(EdgetreeSearch(MHp->togglehead[0], MHp->toggletail[0], nwp[0].outedges) != 0);

      MHp->togglehead[1]=MHp->togglehead[0];
    MHp->toggletail[1]=MHp->toggletail[0];
    MHp->togglehead[0]=head;
    MHp->toggletail[0]=tail;
    nde = nddyads / 2;
    ndn = nddyads / 2;
    nce = nedges-nde;
    ncn = ndyads-nedges-ndn;
/*    MHp->ratio = (nddyads*nddyads) / (odds*(nnodes-nddyads-2)*(nnodes-nddyads-2)); */
    MHp->ratio = (nde*ndn*1.0) / (odds*(nce+1)*(ncn+1));
/*    MHp->ratio = (1.0*(nce-1)*(ncn-1)) / (nde*ndn*odds); */
/*    MHp->ratio = 1.0; */
/*   Rprintf("disconcord nde %d nce %d ndn %d ncn %d nddyads %d MHp->ratio %f\n", */
/*	    nde, nce, ndn,ncn,nddyads, MHp->ratio); */
  }else{
    /* First, select concordant edge at random */
    do{
      rane = 1 + unif_rand() * nedges;
      FindithEdge(MHp->togglehead, MHp->toggletail, rane, &nwp[0]);
    }while(EdgetreeSearch(MHp->togglehead[0], MHp->toggletail[0], nwp[1].outedges) == 0);
       
    /* Next, select concord non-edge at random */
    do{
      head = 1 + unif_rand() * nb1;
      tail = 1 + nb1 + unif_rand() * (nnodes - nb1);
    }while((EdgetreeSearch(head,tail,nwp[0].outedges)!=0) ||
	     (EdgetreeSearch(head,tail,nwp[1].outedges)!=0));

    MHp->togglehead[1]=head;
    MHp->toggletail[1]=tail;
    nde = nddyads / 2;
    ndn = nddyads / 2;
    nce = nedges-nde;
    ncn = ndyads-nedges-ndn;
/*    MHp->ratio = ((nnodes-nddyads)*(nnodes-nddyads)) / (odds*(nddyads+2)*(nddyads+2)); */
    if(nddyads > 4){
      MHp->ratio = (odds*nce*ncn) / ((nde+1)*(ndn+1)*1.0);
/*    MHp->ratio = ((nde+1)*(ndn+1)*odds) / (1.0*nce*ncn); */
    }else{
      MHp->ratio = 100000000.0;
    }
/*   Rprintf("concord nde %d nce %d ndn %d ncn %d nddyads %d MHp->ratio %f\n", */
/*	    nde, nce, ndn,ncn,nddyads, MHp->ratio); */
  }
/*   Rprintf("h0 %d t0 %d h1 %d t1 %d\n", MHp->togglehead[0],  MHp->toggletail[0],  */
/*                                        MHp->togglehead[1],  MHp->toggletail[1]);  */
}

/********************
   void MH_BipartiteHammingTNT
   Tie/no tie:  Gives at least 50% chance of
   proposing a toggle of an existing edge, as opposed
   to simple random toggles that rarely do so in sparse 
   networks
   MSH: The name Hamming is a hack for the Hamming proposals
        It is no different the MH_BipartiteTNT
***********************/
void MH_BipartiteHammingTNT (MHproposal *MHp, DegreeBound *bd, Network *nwp) 
{  
  Vertex head, tail;
  Edge rane, nddyads=nwp[1].nedges;
  int nd, nc;
  static double comp=0.5;
  static double odds;
  static Edge ndyads;
  static Edge nnodes;
  static Edge nb1;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    odds = comp/(1.0-comp);
    nnodes = nwp[0].nnodes;
    nb1 = nwp[0].bipartite;
    ndyads = (nnodes-nb1)*nb1;  
    return;
  }
  
  if (unif_rand() < comp && nddyads > 0) { /* Select a discordant dyad at random */
    rane = 1 + unif_rand() * nddyads;
    FindithEdge(MHp->togglehead, MHp->toggletail, rane, &nwp[1]);
    nd = nddyads;
    nc = ndyads-nd;
    MHp->ratio = (nd*1.0) / (odds*(nc+1));
/*    MHp->ratio = (1.0*(nce-1)*(ncn-1)) / (nde*ndn*odds); */
/*    MHp->ratio = 1.0; */
/*   Rprintf("disconcord nd %d nc %d nddyads %d MHp->ratio %f\n", */
/*	    nd, nc, nddyads, MHp->ratio); */
  }else{
    /* select a concordant dyad at random */
    do{
      head = 1 + unif_rand() * nb1;
      tail = 1 + nb1 + unif_rand() * (nnodes - nb1);
    }while(EdgetreeSearch(head,tail,nwp[1].outedges)!=0);

    MHp->togglehead[0]=head;
    MHp->toggletail[0]=tail;
    nd = nddyads;
    nc = ndyads-nd;
    MHp->ratio = (odds*nc) / ((nd+1)*1.0);
/*   Rprintf("concord nd %d nc %d nddyads %d MHp->ratio %f\n", */
/*	    nd, nc, nddyads, MHp->ratio); */
  }
/*   Rprintf("h0 %d t0 %d h1 %d t1 %d\n", MHp->togglehead[0],  MHp->toggletail[0],  */
/*                                        MHp->togglehead[1],  MHp->toggletail[1]);  */
}


/*********************
 void MH_BipartiteCondDegreeDist
 Pick three nodes -- head, tail, alter -- such that
   * H shares an edge with T
   * H does not share an edge with A
   * deg(T) = deg(A) + 1
 Then, propose swapping the (H,T) edge for a (H,A) edge so that
 deg(H) stays the same while deg(T) and deg(A) swap
 with one another.
*********************/
void MH_BipartiteCondDegreeDist (MHproposal *MHp, DegreeBound *bd, Network *nwp) {  
  int valid, count;
  Edge e;
  Vertex H, T, A, Hin, Hout, Tdeg, Adeg, minA, maxA, i, k;
  double u;
  TreeNode *H_edges;
  Vertex *TA_degrees; 

  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=2;    
    return;
  }
  
  MHp->ratio = 1.0; /* By symmetry:  P(choosing H,T,A) must equal  */
  /*                   P(choosing H,A,T after T and A swap roles), */
  /*                   which makes the ratio equal to 1.           */
  
  for(valid = count = 0; valid == 0 && count<MAX_TRIES/10; count++) {
    Hin = Hout = 0;  
    /* choose a node at random; ensure it has some edges */
    while (Hin + Hout == 0) {
      u = unif_rand();
      if (u < .5) { /* Pick "male" and "female" nodes with equal prob */
        H = 1 + unif_rand() * nwp->bipartite;
      } else {
        H = 1 + nwp->bipartite + unif_rand() * (nwp->nnodes - nwp->bipartite);
      }
      Hin = nwp->indegree[H];
      Hout = nwp->outdegree[H];
    }
    
    /* select an edge to/from H at random */
    k = (int)(unif_rand() * (Hout + Hin));  
    if (k < Hout) { /* we chose an outedge */
      H_edges = nwp->outedges;
      i = k;
      TA_degrees = nwp->indegree;
    } else { /* we chose an inedge */
      H_edges = nwp->inedges;
      i = k - Hout;
      TA_degrees = nwp->outdegree;
    }
    /* Find ith edge in correct edgetree for H; this will be (H,T) or (T,H) */
    e=EdgetreeMinimum(H_edges, H);
    while (i-- > 0) {
      e=EdgetreeSuccessor(H_edges, e);
    } 
    T = H_edges[e].value; 
    Tdeg = nwp->directed_flag ? TA_degrees[T] : nwp->indegree[T] + nwp->outdegree[T];
    
    /* Now choose alter at random */
    /* Now search for eligible alters */
    minA = (u<.5) ? 1 + nwp->bipartite : 1 ;
    maxA = (u<.5) ? nwp->nnodes : nwp->bipartite;    
    i=0;
    for (A=minA; A<=maxA; A++) {
      Adeg = nwp->directed_flag ? TA_degrees[A] : nwp->indegree[A] + nwp->outdegree[A];
      /* To be a valid choice, alter (A) must not be tied to H and it must have */
      /* a degree one less than the tail (T) */
      if (Tdeg == Adeg + 1 && EdgetreeSearch(H, A, H_edges) == 0) {
        if (nwp->directed_flag || EdgetreeSearch(A, H, H_edges) ==0) {
          i++;
        }
      }
    } /* After for-loop, i is # of eligible alters.  */
    if (i>0) {
      valid = 1;
      i = 1 + unif_rand()*i; /* Pick an eligible alter at random */
      for (A=minA; i>0; A++) {
        Adeg = nwp->directed_flag ? TA_degrees[A] : nwp->indegree[A] + nwp->outdegree[A];
        /* To be a valid choice, alter (A) must not be tied to H and it must have */
        /* a degree one less than the tail (T) */
        if (Tdeg == Adeg + 1 && EdgetreeSearch(H, A, H_edges) == 0) {
          if (nwp->directed_flag || EdgetreeSearch(A, H, H_edges) ==0) {
            i--; /* By counting down, when i==0 we have the selected A. */
          }
        }
      } 
    }    
  }
      
  if ( (!nwp->directed_flag && H > T) ||
    (nwp->directed_flag && k < Hout) )
  {
    MHp->togglehead[0] = T;
    MHp->toggletail[0] = H;
  }else{
    MHp->togglehead[0] = H;
    MHp->toggletail[0] = T;
  }

  if(!valid) {
    MHp->togglehead[1] = MHp->togglehead[0];
    MHp->toggletail[1] = MHp->togglehead[0];
  } else {
    if ( (!nwp->directed_flag && H > A) ||
      (nwp->directed_flag && k < Hout) )
    {
      MHp->togglehead[0] = A;
      MHp->toggletail[0] = H;
    }else{
      MHp->togglehead[0] = H;
      MHp->toggletail[0] = A;
    }
  }
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

void MH_BipartiterandomtoggleNonObserved (MHproposal *MHp, DegreeBound *bd, Network *nwp)  {  
  Edge rane;

  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    return;
  }
  MHp->ratio = 1.0;
  rane = 1 + unif_rand() * nwp[1].nedges;
  FindithEdge(MHp->togglehead, MHp->toggletail, rane, &nwp[1]);
/* Rprintf("bip %d nedges %d h %d t %d\n", nwp[1].bipartite, nwp[1].nedges, MHp->togglehead[0],  MHp->toggletail[0]); */
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
