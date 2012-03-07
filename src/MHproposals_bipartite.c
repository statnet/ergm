/*
 *  File ergm/src/MHproposals_bipartite.c
 *  Part of the statnet package, http://statnet.org
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) in
 *    http://statnet.org/attribution
 *
 *  Copyright 2012 the statnet development team
 */
#include "MHproposals_bipartite.h" 

/* Shorthand. */
#define Mtail (MHp->toggletail)
#define Mhead (MHp->togglehead)


/*********************
 void MH_bipartite

 Default MH algorithm for bipartite networks
*********************/
void MH_Bipartiterandomtoggle (MHproposal *MHp, Network *nwp)  {  
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    return;
  }

  /* *** don't forget, edges are (tail, head) now */

  Mtail[0] = 1 + unif_rand() * nwp->bipartite;
  Mhead[0] = 1 + nwp->bipartite + 
                       unif_rand() * (nwp->nnodes - nwp->bipartite);
}

/********************
   void MH_BipartiteConstantEdges
   Chooses a pair of toggles - one a tie and one not. 
***********************/
void MH_BipartiteConstantEdges (MHproposal *MHp, Network *nwp)  {  
  Vertex tail, head;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=2;    
    return;
  } /* Note:  This proposal cannot be used for full or empty observed graphs.
       If desired, we could check for this at initialization phase. 
       (For now, however, no way to easily return an error message and stop.)*/
  /* First, select edge at random */
  GetRandEdge(Mtail, Mhead, nwp);
  /* Second, select dyad at random until it has no edge */

  do{
/*    tail = 1 + unif_rand() * nwp->nnodes; */
/*    head = 1 + unif_rand() * nwp->nnodes; */

  /* *** don't forget, edges are (tail, head) now */

    tail = 1 + unif_rand() * nwp->bipartite;
    head = 1 + nwp->bipartite + 
      unif_rand() * (nwp->nnodes - nwp->bipartite);
  }while(EdgetreeSearch(tail, head, nwp->outedges) != 0);
      
  Mtail[1]=tail;
  Mhead[1]=head;
}

/********************
   void MH_BipartiteTNT
   Tie/no tie:  Gives at least 50% chance of
   proposing a toggle of an existing edge, as opposed
   to simple random toggles that rarely do so in sparse 
   networks
***********************/
void MH_BipartiteTNT (MHproposal *MHp, Network *nwp) 
{
  /* *** don't forget, edges are (tail, head) now */
  Vertex tail, head;
  Edge nedges=nwp->nedges;
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
      GetRandEdge(Mtail, Mhead, nwp);
      /* Thanks to Robert Goudie for pointing out an error in the previous 
      version of this sampler when proposing to go from nedges==0 to nedges==1 
      or vice versa.  Note that this happens extremely rarely unless the 
      network is small or the parameter values lead to extremely sparse 
      networks.  */
      MHp->logratio += log(nedges==1 ? 1.0/(comp*ndyads + (1.0-comp)) :
			 nedges / (odds*ndyads + nedges));
    }else{ /* Select a (tail, head) dyad at random */
      tail = 1 + unif_rand() * nb1;
      head = 1 + nb1 + unif_rand() * (nnodes - nb1);
      Mtail[0] = tail;
      Mhead[0] = head;
      if(EdgetreeSearch(Mtail[0],Mhead[0],nwp->outedges)!=0){
        MHp->logratio += log(nedges==1 ? 1.0/(comp*ndyads + (1.0-comp)) :
			  nedges / (odds*ndyads + nedges));
      }else{
        MHp->logratio += log(nedges==0 ? comp*ndyads + (1.0-comp) :
			  1.0 + (odds*ndyads)/(nedges + 1));
      }
    }
    if(CheckTogglesValid(MHp, nwp)) break;
  }
  
  /* If tries ran out, return failure code. */
  if(trytoggle >= MAX_TRIES){
    Mtail[0]=MH_FAILED;
    Mhead[0]=MH_UNSUCCESSFUL; 
    return;
  }
}

/********************
   void MH_BipartiteHammingConstantEdges
   Chooses a pair of toggles - one a tie and one not. 
   MSH: The name Hamming is a hack for the Hamming proposals
        It is no different the MH_BipartiteConstantEdges
***********************/
void MH_BipartiteHammingConstantEdges (MHproposal *MHp, Network *nwp) 
{
  /* *** don't forget, edges are (tail, head) now */

  Vertex tail, head;
  Edge nedges=nwp[0].nedges, nddyads=nwp[1].nedges;
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
      GetRandEdge(Mtail, Mhead, &nwp[1]);
    }while(EdgetreeSearch(Mtail[0], Mhead[0], nwp[0].outedges) == 0);

    tail=Mtail[0];
    head=Mhead[0];
    /* Next, select discord non-edge at random */
  /* *** don't forget, edges are (tail, head) now */

    do{
      GetRandEdge(Mtail, Mhead, &nwp[1]);
    }while(EdgetreeSearch(Mtail[0], Mhead[0], nwp[0].outedges) != 0);

    Mtail[1]=Mtail[0];
    Mhead[1]=Mhead[0];
    Mtail[0]=tail;
    Mhead[0]=head;
    nde = nddyads / 2;
    ndn = nddyads / 2;
    nce = nedges-nde;
    ncn = ndyads-nedges-ndn;
/*    MHp->ratio = (nddyads*nddyads) / (odds*(nnodes-nddyads-2)*(nnodes-nddyads-2)); */
    MHp->logratio += log((nde*ndn*1.0) / (odds*(nce+1)*(ncn+1)));
/*    MHp->ratio = (1.0*(nce-1)*(ncn-1)) / (nde*ndn*odds); */
/*    MHp->ratio = 1.0; */
/*   Rprintf("disconcord nde %d nce %d ndn %d ncn %d nddyads %d MHp->ratio %f\n", */
/*	    nde, nce, ndn,ncn,nddyads, MHp->ratio); */
  }else{
    /* First, select concordant edge at random */
    do{
      GetRandEdge(Mtail, Mhead, &nwp[0]);
    }while(EdgetreeSearch(Mtail[0], Mhead[0], nwp[1].outedges) == 0);
       
    /* Next, select concord non-edge at random */
    do{
      tail = 1 + unif_rand() * nb1;
      head = 1 + nb1 + unif_rand() * (nnodes - nb1);
    }while((EdgetreeSearch(tail,head,nwp[0].outedges)!=0) ||
	     (EdgetreeSearch(tail,head,nwp[1].outedges)!=0));

    Mtail[1]=tail;
    Mhead[1]=head;
    nde = nddyads / 2;
    ndn = nddyads / 2;
    nce = nedges-nde;
    ncn = ndyads-nedges-ndn;
/*    MHp->ratio = ((nnodes-nddyads)*(nnodes-nddyads)) / (odds*(nddyads+2)*(nddyads+2)); */
    if(nddyads > 4){
      MHp->logratio += log((odds*nce*ncn) / ((nde+1)*(ndn+1)*1.0));
/*    MHp->ratio = ((nde+1)*(ndn+1)*odds) / (1.0*nce*ncn); */
    }else{
      MHp->logratio += 100000000.0;
    }
/*   Rprintf("concord nde %d nce %d ndn %d ncn %d nddyads %d MHp->ratio %f\n", */
/*	    nde, nce, ndn,ncn,nddyads, MHp->ratio); */
  }
/*   Rprintf("h0 %d t0 %d h1 %d t1 %d\n", Mtail[0],  Mhead[0],  */
/*                                        Mtail[1],  Mhead[1]);  */
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
void MH_BipartiteHammingTNT (MHproposal *MHp, Network *nwp) 
{
  /* *** don't forget, edges are (tail, head) now */  
  Vertex tail, head;
  Edge nddyads=nwp[1].nedges;
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
    GetRandEdge(Mtail, Mhead, &nwp[1]);
    nd = nddyads;
    nc = ndyads-nd;
    /*  Fixme!  Not sure whether the ratio is calculated correctly here.
        Check out the similar ratio calculations for other TNT proposals. */
    MHp->logratio += log((nd*1.0) / (odds*(nc+1)));
/*    MHp->ratio = (1.0*(nce-1)*(ncn-1)) / (nde*ndn*odds); */
/*    MHp->ratio = 1.0; */
/*   Rprintf("disconcord nd %d nc %d nddyads %d MHp->ratio %f\n", */
/*	    nd, nc, nddyads, MHp->ratio); */
  }else{
    /* select a concordant dyad at random */
    do{
      tail = 1 + unif_rand() * nb1;
      head = 1 + nb1 + unif_rand() * (nnodes - nb1);
    }while(EdgetreeSearch(tail,head,nwp[1].outedges)!=0);

    Mtail[0]=tail;
    Mhead[0]=head;
    nd = nddyads;
    nc = ndyads-nd;
    /*  Fixme!  Not sure whether the ratio is calculated correctly here.
        Check out the similar ratio calculations for other TNT proposals. */
    MHp->logratio += log((odds*nc) / ((nd+1)*1.0));
/*   Rprintf("concord nd %d nc %d nddyads %d MHp->ratio %f\n", */
/*	    nd, nc, nddyads, MHp->ratio); */
  }
/*   Rprintf("h0 %d t0 %d h1 %d t1 %d\n", Mtail[0],  Mhead[0],  */
/*                                        Mtail[1],  Mhead[1]);  */
}


/*********************
 void MH_BipartiteCondDegreeDist
 Pick three nodes -- tail, head, alter -- such that
   * tail shares an edge with head
   * tail does not share an edge with A
   * deg(head) = deg(A) + 1
 Then, propose swapping the (tail,head) edge for a (tail,A) edge so that
 deg(tail) stays the same while deg(head) and deg(A) swap
 with one another.
*********************/
void MH_BipartiteCondDegreeDist (MHproposal *MHp, Network *nwp) {  
  int valid, count;
  Edge e;
  Vertex tail, head, A, tailin, tailout, headdeg, Adeg, minA, maxA, i, k;
  double u;
  TreeNode *tail_edges;
  Vertex *headA_degrees; 

  /* *** don't forget, edges are (tail, head) now */

  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=2;    
    return;
  }
  
  MHp->logratio += 0;  /* By symmetry:  P(choosing tail,head,A) must equal  */
  /*                   P(choosing tail,A,head after head and A swap roles), */
  /*                   which makes the ratio equal to 1.           */
  
  for(valid = count = 0; valid == 0 && count<MAX_TRIES/10; count++) {
    tailin = tailout = 0;  
    /* choose a node at random; ensure it has some edges */
    while (tailin + tailout == 0) {
      u = unif_rand();
      if (u < .5) { /* Pick "male" and "female" nodes with equal prob */
        tail = 1 + unif_rand() * nwp->bipartite;
      } else {
        tail = 1 + nwp->bipartite + unif_rand() * (nwp->nnodes - nwp->bipartite);
      }
      tailin = nwp->indegree[tail];
      tailout = nwp->outdegree[tail];
    }
    
    /* select an edge to/from tail at random */
    k = (int)(unif_rand() * (tailout + tailin));  
    if (k < tailout) { /* we chose an outedge */
      tail_edges = nwp->outedges;
      i = k;
      headA_degrees = nwp->indegree;
    } else { /* we chose an inedge */
      tail_edges = nwp->inedges;
      i = k - tailout;
      headA_degrees = nwp->outdegree;
    }
    /* Find ith edge in correct edgetree for tail; this will be (tail,head) or (head,tail) */
    e=EdgetreeMinimum(tail_edges, tail);
    while (i-- > 0) {
      e=EdgetreeSuccessor(tail_edges, e);
    } 
    head = tail_edges[e].value; 
    headdeg = nwp->directed_flag ? headA_degrees[head] : nwp->indegree[head] + nwp->outdegree[head];
    
    /* Now choose alter at random */
    /* Now search for eligible alters */
    minA = (u<.5) ? 1 + nwp->bipartite : 1 ;
    maxA = (u<.5) ? nwp->nnodes : nwp->bipartite;    
    i=0;
    for (A=minA; A<=maxA; A++) {
      Adeg = nwp->directed_flag ? headA_degrees[A] : nwp->indegree[A] + nwp->outdegree[A];
      /* To be a valid choice, alter (A) must not be tied to tail and it must have */
      /* a degree one less than the head */
      if (headdeg == Adeg + 1 && EdgetreeSearch(tail, A, tail_edges) == 0) {
        if (nwp->directed_flag || EdgetreeSearch(A, tail, tail_edges) ==0) {
          i++;
        }
      }
    } /* After for-loop, i is # of eligible alters.  */
    if (i>0) {
      valid = 1;
      i = 1 + unif_rand()*i; /* Pick an eligible alter at random */
      for (A=minA; i>0; A++) {
        Adeg = nwp->directed_flag ? headA_degrees[A] : nwp->indegree[A] + nwp->outdegree[A];
        /* To be a valid choice, alter (A) must not be tied to tail and it must have */
        /* a degree one less than the head */
        if (headdeg == Adeg + 1 && EdgetreeSearch(tail, A, tail_edges) == 0) {
          if (nwp->directed_flag || EdgetreeSearch(A, tail, tail_edges) ==0) {
            i--; /* By counting down, when i==0 we have the selected A. */
          }
        }
      } 
    }    
  }
      
  if ( (!nwp->directed_flag && tail > head) ||
    (nwp->directed_flag && k < tailout) )
  {
    Mtail[0] = head;
    Mhead[0] = tail;
  }else{
    Mtail[0] = tail;
    Mhead[0] = head;
  }

  if(!valid) {
    Mtail[1] = Mtail[0];
    Mhead[1] = Mtail[0];
  } else {
    if ( (!nwp->directed_flag && tail > A) ||
      (nwp->directed_flag && k < tailout) )
    {
      Mtail[0] = A;
      Mhead[0] = tail;
    }else{
      Mtail[0] = tail;
      Mhead[0] = A;
    }
  }
}


void MH_BipartiterandomtoggleNonObserved (MHproposal *MHp, Network *nwp)  {  
  Edge rane, nmissing = MHp->inputs[0];

  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    return;
  }

  // Note that missing edgelist is indexed from 0 but the first
  // element of MHp->inputs is the number of missing edges.
  rane = 1 + unif_rand() * nmissing;
  
  Mtail[0]=MHp->inputs[rane];
  Mhead[0]=MHp->inputs[nmissing+rane];
}

/* CondDegree */
void MH_BipartiteCondDegTetradToggles (MHproposal *MHp, Network *nwp)  {  
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
  /* First, select edge at random */
  GetRandEdge(&A, &B, nwp);
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
  Mtail[0]=A; Mhead[0]=B;
  Mtail[1]=tmpA; Mhead[1]=tmpC;
  Mtail[2]=tmpD; Mhead[2]=tmpB;
  Mtail[3]=MIN(C,D); Mhead[3]=MAX(C,D);
   Rprintf("tail0 %d head1 %d\n",Mtail[0], Mhead[0]);
   Rprintf("tail0 %d head1 %d\n",Mtail[1], Mhead[1]);
   Rprintf("tail0 %d head1 %d\n",Mtail[2], Mhead[2]);
   Rprintf("tail0 %d head1 %d\n",Mtail[3], Mhead[3]);
}  
  
void MH_BipartiteCondDegree (MHproposal *MHp, Network *nwp)  {  
  
  if(MHp->ntoggles == 0) { /* Initialize CondDeg by */
	                   /* Choosing Hexad or Tetrad */
    MH_BipartiteCondDegHexadToggles (MHp, nwp);
    MH_BipartiteCondDegTetradToggles (MHp, nwp);
    MHp->ntoggles=4;
    return;
  }

  if( unif_rand() > 0.9 ){
    MHp->ntoggles=6;
  }else{
    MHp->ntoggles=4;
  }

  if(MHp->ntoggles == 6) { /* Call Hexad */
    MH_BipartiteCondDegHexadToggles (MHp, nwp);
  }else{ /* call Tetrad */
    MH_BipartiteCondDegTetradToggles (MHp, nwp);
  }
}
void MH_BipartiteCondDegHexadToggles (MHproposal *MHp, Network *nwp)  {  
  int x1, x2, x3, x4, x5, x6;
  int fvalid, trynode;
  Vertex tail1, tail2, tail3, head1, head2, head3;
  static Edge nnodes;
  static Edge nb1, nb2;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=6;    
    nnodes = nwp[0].nnodes;
    nb1 = nwp[0].bipartite;
    nb2 = nnodes-nb1;
    return;
  }
  
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
  
  tail1 = 1 + nb1 + unif_rand() * nb2;
  head2 = 1 + unif_rand() * nb1;
  head2 = 1 + unif_rand() * nb1;
  head3 = 1 + unif_rand() * nb1;
  while(head3 == head2){
    head3 = 1 + unif_rand() * nb1;
  }
  if ((!nwp->directed_flag) && tail1 > head2){
    x1 = EdgetreeSearch(head2, tail1, nwp->outedges) > 0;
  }else{
    x1 = EdgetreeSearch(tail1, head2, nwp->outedges) > 0;
  }
  if ((!nwp->directed_flag) && tail1 > head3){
    x2 = EdgetreeSearch(head3, tail1, nwp->outedges) > 0;
  }else{
    x2 = EdgetreeSearch(tail1, head3, nwp->outedges) > 0;
  }
  if (x1 != x2){
    head1 = 1 + unif_rand() * nb1;
    while(head1 == head2 || head1 == head3 || head1 == tail1){
      head1 = 1 + unif_rand() * nb1;
    }
    tail2 = 1 + nb1 + unif_rand() * nb2;
    while(tail2 == head2 || tail2 == head3 || tail2 == tail1 ||
	  tail2 == head1 ){
      tail2 = 1 + nb1 + unif_rand() * nb2;
    }
    if ((!nwp->directed_flag) && tail2 > head1){
      x3 = EdgetreeSearch(head1, tail2, nwp->outedges) > 0;
    }else{
      x3 = EdgetreeSearch(tail2, head1, nwp->outedges) > 0;
    }
    if (x2 == x3){
      if ((!nwp->directed_flag) && tail2 > head3){
	x4 = EdgetreeSearch(head3, tail2, nwp->outedges) > 0;
      }else{
	x4 = EdgetreeSearch(tail2, head3, nwp->outedges) > 0;
      }
      if (x4 == x1){
	tail3 = 1 + unif_rand() * nb2;
	while(tail3 == head2 || tail3 == head3 || tail3 == tail1 ||
	      tail3 == head1 || tail3 == tail2 ){
	  tail3 = 1 + nb1 + unif_rand() * nb2;
	}
	if ((!nwp->directed_flag) && tail3 > head1){
	  x5 = EdgetreeSearch(head1, tail3, nwp->outedges) > 0;
	}else{
	  x5 = EdgetreeSearch(tail3, head1, nwp->outedges) > 0;
	}
	if (x5 == x1){
	  if ((!nwp->directed_flag) && tail3 > head2){
	    x6 = EdgetreeSearch(head2, tail3, nwp->outedges) > 0;
	  }else{
	    x6 = EdgetreeSearch(tail3, head2, nwp->outedges) > 0;
	  }
	  if (x6 == x2){
	    if ( (!nwp->directed_flag) ){
	      if ( tail1 > head2 ){
		Mtail[0] = head2;
		Mhead[0] = tail1;
	      }else{
		Mtail[0] = tail1;
		Mhead[0] = head2;
	      }
	      if ( tail1 > head3 ){
		Mtail[1] = head3;
		Mhead[1] = tail1;
	      }else{
		Mtail[1] = tail1;
		Mhead[1] = head3;
	      }
	      if ( tail2 > head1 ){
		Mtail[2] = head1;
		Mhead[2] = tail2;
	      }else{
		Mtail[2] = tail2;
		Mhead[2] = head1;
	      }
	      if ( tail2 > head3 ){
		Mtail[3] = head3;
		Mhead[3] = tail2;
	      }else{
		Mtail[3] = tail2;
		Mhead[3] = head3;
	      }
	      if ( tail3 > head1 ){
		Mtail[4] = head1;
		Mhead[4] = tail3;
	      }else{
		Mtail[4] = tail3;
		Mhead[4] = head1;
	      }
	      if ( tail3 > head2 ){
		Mtail[5] = head2;
		Mhead[5] = tail3;
	      }else{
		Mtail[5] = tail3;
		Mhead[5] = head2;
	      }
	    }else{
	      Mtail[0] = tail1;
	      Mhead[0] = head2;
	      Mtail[1] = tail1;
	      Mhead[1] = head3;
	      Mtail[2] = tail2;
	      Mhead[2] = head1;
	      Mtail[3] = tail2;
	      Mhead[3] = head3;
	      Mtail[4] = tail3;
	      Mhead[4] = head1;
	      Mtail[5] = tail3;
	      Mhead[5] = head2;
	    }
	    fvalid = 1;
	  }
	}
      }
    }
  }
  }

  if(trynode==5000){
      Mtail[0] = 1;
      Mhead[0] = 2;
      Mtail[1] = 1;
      Mhead[1] = 2;
      Mtail[2] = 1;
      Mhead[2] = 2;
      Mtail[3] = 1;
      Mhead[3] = 2;
      Mtail[4] = 1;
      Mhead[4] = 2;
      Mtail[5] = 1;
      Mhead[5] = 2;
  }
}
