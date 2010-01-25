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

/* CondDegree */
void MH_BipartiteCondDegTetradToggles (MHproposal *MHp, DegreeBound *bd, Network *nwp)  {  
  Vertex A, B, C, D=0;
  Vertex tmpA=0, tmpB, tmpC, tmpD=0;
  int valid, n_C_nbrs, i;
  TreeNode *tree;
  Edge e;
  static Edge nnodes;
  static Edge nb1, nb2;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=4;    
    nnodes = nwp[0].nnodes;
    nb1 = nwp[0].bipartite;
    nb2 = nnodes-nb1;
   Rprintf("nb1 %d nb2 %d \n",nb1,nb2);
    return;
  } /* Note:  This proposal does not make work (well) with 
  directed graphs; however, we haven't yet implemented a way
  to warn the user about this.  */
  MHp->ratio=1.0;
  /* First, select edge at random */
  FindithEdge(&A, &B, 1+nwp->nedges*unif_rand(), nwp);
   Rprintf("nb1 %d nb2 %d A %d B %d\n",nb1,nb2,A,B);
  /* Second, select a non-neighbor C of B and a random neighbor
  D of C such that D is not a neighbor of A.  */
  valid=0;
  while (!valid) {
    C = 1+nb1+unif_rand() * nb2;
    tmpB = MIN(B,C);
    tmpC = MAX(B,C);
   Rprintf("nb1 %d nb2 %d A %d B %d C %d\n",nb1,nb2,A,B,C);
    if (C != A && C != B  && EdgetreeSearch(tmpB, tmpC, nwp->outedges)==0) {
      /* Now pick D, a random neighbor of C, */
      n_C_nbrs = nwp->indegree[C]+nwp->outdegree[C];
   Rprintf("nb1 %d nb2 %d A %d B %d C %d nC %d idc %d\n",nb1,nb2,A,B,C,n_C_nbrs, nwp->indegree[C]);
      if(n_C_nbrs > 0) {
        i = 1 + unif_rand() * n_C_nbrs;
        tree = nwp->inedges;
        if (i < nwp->indegree[C]) {
          i = i - nwp->indegree[C];
          tree = nwp->inedges;
        }
        e=EdgetreeMinimum(tree, C);
        while (i-- > 1) {
          e=EdgetreeSuccessor(tree,e);
        }
        D = tree[e].value;
   Rprintf("nb1 %d nb2 %d A %d B %d C %d, nC %d D %d\n",nb1,nb2,A,B,C,n_C_nbrs,D);
        /* Now check to ensure that (A,D) does not exist */        
        tmpA = MIN(A,D);
        tmpD = MAX(A,D);
   Rprintf("tmpA %d tmpB %d\n",tmpA, tmpB);
        if (A !=D && B != D && EdgetreeSearch(tmpA, tmpD, nwp->outedges)==0) {
          valid=1;
        }
      }
    }
  }
  MHp->togglehead[0]=A; MHp->toggletail[0]=B;
  MHp->togglehead[1]=tmpA; MHp->toggletail[1]=tmpC;
  MHp->togglehead[2]=tmpD; MHp->toggletail[2]=tmpB;
  MHp->togglehead[3]=MIN(C,D); MHp->toggletail[3]=MAX(C,D);
   Rprintf("h0 %d t1 %d\n",MHp->togglehead[0], MHp->toggletail[0]);
   Rprintf("h0 %d t1 %d\n",MHp->togglehead[1], MHp->toggletail[1]);
   Rprintf("h0 %d t1 %d\n",MHp->togglehead[2], MHp->toggletail[2]);
   Rprintf("h0 %d t1 %d\n",MHp->togglehead[3], MHp->toggletail[3]);
}  
  
void MH_BipartiteCondDegree (MHproposal *MHp, DegreeBound *bd, Network *nwp)  {  
  MHp->ratio=1.0;
  
  if(MHp->ntoggles == 0) { /* Initialize CondDeg by */
	                   /* Choosing Hexad or Tetrad */
    MH_BipartiteCondDegHexadToggles (MHp, bd, nwp);
    MH_BipartiteCondDegTetradToggles (MHp, bd, nwp);
    MHp->ntoggles=4;
    return;
  }

  if( unif_rand() > 0.9 ){
    MHp->ntoggles=6;
  }else{
    MHp->ntoggles=4;
  }

  if(MHp->ntoggles == 6) { /* Call Hexad */
    MH_BipartiteCondDegHexadToggles (MHp, bd, nwp);
  }else{ /* call Tetrad */
    MH_BipartiteCondDegTetradToggles (MHp, bd, nwp);
  }
}
void MH_BipartiteCondDegHexadToggles (MHproposal *MHp, DegreeBound *bd, Network *nwp)  {  
  int x1, x2, x3, x4, x5, x6;
  int fvalid, trynode;
  Vertex head1, head2, head3, tail1, tail2, tail3;
  static Edge nnodes;
  static Edge nb1, nb2;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=6;    
    nnodes = nwp[0].nnodes;
    nb1 = nwp[0].bipartite;
    nb2 = nnodes-nb1;
    return;
  }
  MHp->ratio=1.0;
  
  x1 = -1;
  x2 = -1;
  x3 = -1;
  x4 = -1;
  x5 = -1;
  x6 = -1;

  fvalid = 0;
  trynode = 0;
  while(fvalid==0 && trynode < MAX_TRIES){

  trynode++;
  /* select a node at random */
  
  head1 = 1 + nb1 + unif_rand() * nb2;
  tail2 = 1 + unif_rand() * nb1;
  tail2 = 1 + unif_rand() * nb1;
  tail3 = 1 + unif_rand() * nb1;
  while(tail3 == tail2){
    tail3 = 1 + unif_rand() * nb1;
  }
  if ((!nwp->directed_flag) && head1 > tail2){
    x1 = EdgetreeSearch(tail2, head1, nwp->outedges) > 0;
  }else{
    x1 = EdgetreeSearch(head1, tail2, nwp->outedges) > 0;
  }
  if ((!nwp->directed_flag) && head1 > tail3){
    x2 = EdgetreeSearch(tail3, head1, nwp->outedges) > 0;
  }else{
    x2 = EdgetreeSearch(head1, tail3, nwp->outedges) > 0;
  }
  if (x1 != x2){
    tail1 = 1 + unif_rand() * nb1;
    while(tail1 == tail2 || tail1 == tail3 || tail1 == head1){
      tail1 = 1 + unif_rand() * nb1;
    }
    head2 = 1 + nb1 + unif_rand() * nb2;
    while(head2 == tail2 || head2 == tail3 || head2 == head1 ||
	  head2 == tail1 ){
      head2 = 1 + nb1 + unif_rand() * nb2;
    }
    if ((!nwp->directed_flag) && head2 > tail1){
      x3 = EdgetreeSearch(tail1, head2, nwp->outedges) > 0;
    }else{
      x3 = EdgetreeSearch(head2, tail1, nwp->outedges) > 0;
    }
    if (x2 == x3){
      if ((!nwp->directed_flag) && head2 > tail3){
	x4 = EdgetreeSearch(tail3, head2, nwp->outedges) > 0;
      }else{
	x4 = EdgetreeSearch(head2, tail3, nwp->outedges) > 0;
      }
      if (x4 == x1){
	head3 = 1 + unif_rand() * nb2;
	while(head3 == tail2 || head3 == tail3 || head3 == head1 ||
	      head3 == tail1 || head3 == head2 ){
	  head3 = 1 + nb1 + unif_rand() * nb2;
	}
	if ((!nwp->directed_flag) && head3 > tail1){
	  x5 = EdgetreeSearch(tail1, head3, nwp->outedges) > 0;
	}else{
	  x5 = EdgetreeSearch(head3, tail1, nwp->outedges) > 0;
	}
	if (x5 == x1){
	  if ((!nwp->directed_flag) && head3 > tail2){
	    x6 = EdgetreeSearch(tail2, head3, nwp->outedges) > 0;
	  }else{
	    x6 = EdgetreeSearch(head3, tail2, nwp->outedges) > 0;
	  }
	  if (x6 == x2){
	    if ( (!nwp->directed_flag) ){
	      if ( head1 > tail2 ){
		MHp->togglehead[0] = tail2;
		MHp->toggletail[0] = head1;
	      }else{
		MHp->togglehead[0] = head1;
		MHp->toggletail[0] = tail2;
	      }
	      if ( head1 > tail3 ){
		MHp->togglehead[1] = tail3;
		MHp->toggletail[1] = head1;
	      }else{
		MHp->togglehead[1] = head1;
		MHp->toggletail[1] = tail3;
	      }
	      if ( head2 > tail1 ){
		MHp->togglehead[2] = tail1;
		MHp->toggletail[2] = head2;
	      }else{
		MHp->togglehead[2] = head2;
		MHp->toggletail[2] = tail1;
	      }
	      if ( head2 > tail3 ){
		MHp->togglehead[3] = tail3;
		MHp->toggletail[3] = head2;
	      }else{
		MHp->togglehead[3] = head2;
		MHp->toggletail[3] = tail3;
	      }
	      if ( head3 > tail1 ){
		MHp->togglehead[4] = tail1;
		MHp->toggletail[4] = head3;
	      }else{
		MHp->togglehead[4] = head3;
		MHp->toggletail[4] = tail1;
	      }
	      if ( head3 > tail2 ){
		MHp->togglehead[5] = tail2;
		MHp->toggletail[5] = head3;
	      }else{
		MHp->togglehead[5] = head3;
		MHp->toggletail[5] = tail2;
	      }
	    }else{
	      MHp->togglehead[0] = head1;
	      MHp->toggletail[0] = tail2;
	      MHp->togglehead[1] = head1;
	      MHp->toggletail[1] = tail3;
	      MHp->togglehead[2] = head2;
	      MHp->toggletail[2] = tail1;
	      MHp->togglehead[3] = head2;
	      MHp->toggletail[3] = tail3;
	      MHp->togglehead[4] = head3;
	      MHp->toggletail[4] = tail1;
	      MHp->togglehead[5] = head3;
	      MHp->toggletail[5] = tail2;
	    }
	    fvalid = 1;
	  }
	}
      }
    }
  }
  }

  if(trynode==5000){
      MHp->togglehead[0] = 1;
      MHp->toggletail[0] = 2;
      MHp->togglehead[1] = 1;
      MHp->toggletail[1] = 2;
      MHp->togglehead[2] = 1;
      MHp->toggletail[2] = 2;
      MHp->togglehead[3] = 1;
      MHp->toggletail[3] = 2;
      MHp->togglehead[4] = 1;
      MHp->toggletail[4] = 2;
      MHp->togglehead[5] = 1;
      MHp->toggletail[5] = 2;
  }
}
