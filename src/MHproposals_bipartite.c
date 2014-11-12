/*  File src/MHproposals_bipartite.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2014 Statnet Commons
 */
#include "MHproposals_bipartite.h" 

/* Shorthand. */
#define Mtail (MHp->toggletail)
#define Mhead (MHp->togglehead)


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
  static Dyad ndyads;
  static Edge nnodes;
  static Edge nb1;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=2;
    odds = comp/(1.0-comp);
    nnodes = nwp[0].nnodes;
    nb1 = nwp[0].bipartite;
    ndyads = DYADCOUNT(nnodes, nb1, 0);
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
  static Dyad ndyads;
  static Edge nnodes;
  static Edge nb1;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    odds = comp/(1.0-comp);
    nnodes = nwp[0].nnodes;
    nb1 = nwp[0].bipartite;
    ndyads = DYADCOUNT(nnodes, nb1, 0);
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

