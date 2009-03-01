/*
 *  File ergm/src/MHproposals.c
 *  Part of the statnet package, http://statnetproject.org
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) in
 *    http://statnetproject.org/attribution
 *
 * Copyright 2003 Mark S. Handcock, University of Washington
 *                David R. Hunter, Penn State University
 *                Carter T. Butts, University of California - Irvine
 *                Steven M. Goodreau, University of Washington
 *                Martina Morris, University of Washington
 * Copyright 2007 The statnet Development Team
 */
#include "MHproposals.h"

/* Shorthand. */
#define Mhead (MHp->togglehead)
#define Mtail (MHp->toggletail)

/*********************
 void MH_randomtoggle

 Default MH algorithm
*********************/
void MH_randomtoggle (MHproposal *MHp, DegreeBound *bd, Network *nwp)  {  
  Vertex head, tail;
  int fvalid, trytoggle;
  
  if(MHp->ntoggles == 0) { /* Initialize randomtoggle */
    MHp->ntoggles=1;
    return;
  }
  MHp->ratio = 1.0;
  
  fvalid = 0;
  trytoggle = 0;
  while(fvalid==0 && trytoggle < MAX_TRIES){
    head = 1 + unif_rand() * nwp->nnodes;
    while ((tail = 1 + unif_rand() * nwp->nnodes) == head);
    if (!nwp->directed_flag && head > tail) {
      MHp->togglehead[0] = tail;
      MHp->toggletail[0] = head;
    }else{
      MHp->togglehead[0] = head;
      MHp->toggletail[0] = tail;
    }
    fvalid=CheckTogglesValid(MHp, bd, nwp);
  }
}

/********************
   void MH_TNT
   Tie/no tie:  Gives at least 50% chance of
   proposing a toggle of an existing edge, as opposed
   to simple random toggles that rarely do so in sparse 
   networks
***********************/
void MH_TNT (MHproposal *MHp, DegreeBound *bd, Network *nwp) 
{  
  Vertex head, tail;
  Edge rane, nedges=nwp->nedges;
  static double comp=0.5;
  static double odds;
  static Edge ndyads;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    odds = comp/(1.0-comp);
    ndyads = (nwp->nnodes-1)*nwp->nnodes / (nwp->directed_flag? 1:2);  
    return;
  }
  
  for(int trytoggle = 0; trytoggle < MAX_TRIES; trytoggle++){
    if (unif_rand() < comp && nedges > 0) { /* Select a tie at random */
      rane = 1 + unif_rand() * nedges;
      FindithEdge(MHp->togglehead, MHp->toggletail, rane, nwp);
      MHp->ratio = nedges  / (odds*ndyads + nedges);
    }else{ /* Select a dyad at random */
      do{
        head = 1 + unif_rand() * nwp->nnodes;
      }while ((tail = 1 + unif_rand() * nwp->nnodes) == head);
      if (head > tail && !nwp->directed_flag)  {
        MHp->togglehead[0] = tail;
        MHp->toggletail[0] = head;
      }else{
        MHp->togglehead[0] = head;
        MHp->toggletail[0] = tail;
      }
      if(EdgetreeSearch(MHp->togglehead[0],MHp->toggletail[0],nwp->outedges)!=0){
        MHp->ratio = nedges / (odds*ndyads + nedges);
      }else{
        MHp->ratio = 1.0 + (odds*ndyads)/(nedges + 1);
      }
    }
    if(CheckTogglesValid(MHp,bd,nwp)) break;
  }
}

void MH_TNT10 (MHproposal *MHp, DegreeBound *bd, Network *nwp) 
{  
  Vertex head, tail;
  Edge rane, nedges=nwp->nedges;
  static double comp=0.5;
  static double odds;
  static Edge ndyads;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=10;
    MHp->ratio = 1.0;
    odds = comp/(1.0-comp);
    ndyads = (nwp->nnodes-1)*nwp->nnodes / (nwp->directed_flag? 1:2);  
    return;
  }
  
  for(int trytoggle = 0; trytoggle < MAX_TRIES; trytoggle++){
   for(int n = 0; n < 10; n++){
    if (unif_rand() < comp && nedges > 0) { /* Select a tie at random */
      rane = 1 + unif_rand() * nedges;
      FindithEdge(MHp->togglehead, MHp->toggletail, rane, nwp);
      MHp->ratio *= nedges  / (odds*ndyads + nedges);
    }else{ /* Select a dyad at random */
      do{
        head = 1 + unif_rand() * nwp->nnodes;
      }while ((tail = 1 + unif_rand() * nwp->nnodes) == head);
      if (head > tail && !nwp->directed_flag)  {
        MHp->togglehead[n] = tail;
        MHp->toggletail[n] = head;
      }else{
        MHp->togglehead[n] = head;
        MHp->toggletail[n] = tail;
      }
      if(EdgetreeSearch(MHp->togglehead[n],MHp->toggletail[n],nwp->outedges)!=0){
        MHp->ratio *= nedges / (odds*ndyads + nedges);
      }else{
        MHp->ratio *= 1.0 + (odds*ndyads)/(nedges + 1);
      }
    } 
   }
   if(CheckTogglesValid(MHp,bd,nwp)) break;
  }
}

/*********************
 void MH_constantedges
 propose pairs of toggles that keep number of edges
 the same.  This is done by (a) choosing an existing edge
 at random; (b) repeatedly choosing dyads at random until 
 one is found that does not have an edge; and (c) proposing
 toggling both these dyads.  Note that step (b) will be very 
 inefficient if the network is nearly complete, so this proposal is
 NOT recommended for such networks.  However, most network
 datasets are sparse, so this is not likely to be an issue.
*********************/
void MH_ConstantEdges (MHproposal *MHp, DegreeBound *bd, Network *nwp)  {  
  Vertex head, tail, temp;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=2;    
    MHp->ratio=1.0;   
    return;
  } /* Note:  This proposal cannot be used for full or empty observed graphs.
       If desired, we could check for this at initialization phase. 
       (For now, however, no way to easily return an error message and stop.)*/
  for(int trytoggle = 0; trytoggle < MAX_TRIES; trytoggle++){
    /* First, select edge at random */
    FindithEdge(MHp->togglehead, MHp->toggletail, 1+nwp->nedges*unif_rand(), nwp);
    /* Second, select dyad at random until it has no edge */
    do{
      head = 1 + unif_rand() * nwp->nnodes;
      tail = 1 + unif_rand() * nwp->nnodes;
      if (!nwp->directed_flag && head > tail) {
        temp=head;  head=tail;  tail=temp; /* swap head for tail */
      }
    }while (EdgetreeSearch(head, tail, nwp->outedges) != 0 || head == tail);
    
    MHp->togglehead[1]=head;
    MHp->toggletail[1]=tail;
    
    if(CheckTogglesValid(MHp,bd,nwp)) break; 
  }
}
  
/*********************
 void MH_CondDegTetra
 
 Select an edge (A,B) at random.  Then select
 another node C at random, making sure it has at least one neighbor
 and that C is not a neighbor of B.  Select one of the neighbors of C,
 D, at random.  Thus, we have A--B, C--D, but not B--C.  Finally,
 check to see whether A--D.  If not, then change
      A--B            A  B
                to     \/
                       /\
      C--D            C  D
  
 Note that this algorithm may be inefficient if the network is not sparse. 
*********************/
void MH_CondDegTetra (MHproposal *MHp, DegreeBound *bd, Network *nwp)  {  
  Vertex A, B, C, D=0;
  Vertex tmpA=0, tmpB, tmpC, tmpD=0;
  int valid, n_C_nbrs, i;
  TreeNode *tree;
  Edge e;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=4;    
    return;
  } /* Note:  This proposal does not make work (well) with 
  directed graphs; however, we haven't yet implemented a way
  to warn the user about this.  */
  MHp->ratio=1.0;
  /* First, select edge at random */
  FindithEdge(&A, &B, 1+nwp->nedges*unif_rand(), nwp);
  /* Second, select a non-neighbor C of B and a random neighbor
  D of C such that D is not a neighbor of A.  */
  valid=0;
  while (!valid) {
    C = 1+unif_rand() * nwp->nnodes;
    tmpB = MIN(B,C);
    tmpC = MAX(B,C);
    if (C != A && C != B  && EdgetreeSearch(tmpB, tmpC, nwp->outedges)==0) {
      /* Now pick D, a random neighbor of C, */
      n_C_nbrs = nwp->indegree[C]+nwp->outdegree[C];
      if(n_C_nbrs > 0) {
        i = 1 + unif_rand() * n_C_nbrs;
        tree = nwp->inedges;
        if (i < nwp->indegree[C]) {
          i = i - nwp->indegree[C];
          tree = nwp->outedges;
        }
        e=EdgetreeMinimum(tree, C);
        while (i-- > 1) {
          e=EdgetreeSuccessor(tree,e);
        }
        D = tree[e].value;
        /* Now check to ensure that (A,D) does not exist */        
        tmpA = MIN(A,D);
        tmpD = MAX(A,D);
        if (A !=D && B != D && EdgetreeSearch(tmpA, tmpD, nwp->outedges)==0) {
          valid=1;
        }
      }
    }
  }
  MHp->togglehead[0]=A; MHp->toggletail[0]=B;
  MHp->togglehead[1]=tmpA; MHp->toggletail[1]=tmpD;
  MHp->togglehead[2]=tmpB; MHp->toggletail[2]=tmpC;
  MHp->togglehead[3]=MIN(C,D); MHp->toggletail[1]=MAX(C,D);
}  
/*  I still have some questions about the MH_conddeg routine:
  
  1.  Should we write this so it at could "work" with digraphs?
  2.  Should we have a special function to find the ith edge 
      incident to a given node?
  3.  Does this algorithm introduce any selection bias that should
      be corrected by MHp->ratio?  In other words, is it really true
      that this algorithm gives the same chance for the reverse change
      as for the forward change?
*/

/*********************
 void MH_CondDegreeDist
 It used to be called  MH_CondDegDistSwapToggles
*********************/
void MH_CondDegreeDist (MHproposal *MHp, DegreeBound *bd, Network *nwp) {  
  int noutedge=0, ninedge=0, k, fvalid;
  int k0, j0, j1, k1;
  int j0h, j1h;
  int trynode;
  Vertex e, alter, head=0, tail, tail1;
  MHp->ratio=1.0;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=2;    
    return;
  }

  fvalid = 0;
  trynode = 0;
  while(fvalid==0 && trynode < 500){

  trynode++;
  /* select a node at random */
  while(noutedge+ninedge==0){
    /* select a node at random */
    head = 1 + unif_rand() * nwp->nnodes;
    ninedge  = nwp->indegree[head];
    noutedge = nwp->outdegree[head];
  }

  /* choose a edge of the node at random */
  
  k0 = (int)(unif_rand() * (noutedge+ninedge)); 
  if (k0 < noutedge){
    k=0;
    for(e = EdgetreeMinimum(nwp->outedges, head);
    ((tail = nwp->outedges[e].value) != 0 && k<k0);
    e = EdgetreeSuccessor(nwp->outedges, e)){++k;}
  }else{
    k=0;
    for(e = EdgetreeMinimum(nwp->inedges, head);
    ((tail = nwp->inedges[e].value) != 0 && k<(k0-noutedge));
    e = EdgetreeSuccessor(nwp->inedges, e)){++k;}
  }

  if ( (!nwp->directed_flag && head > tail) ||
  (nwp->directed_flag && k0 >= noutedge) ) {
    MHp->togglehead[0] = tail;
    MHp->toggletail[0] = head;
  }else{
    MHp->togglehead[0] = head;
    MHp->toggletail[0] = tail;
  }
  
  k1=0;
  fvalid=0;
  while(fvalid==0 && k1 < 100){
    while((alter = 1 + unif_rand() * nwp->nnodes) == head);
    fvalid=1;
    if(alter == tail){fvalid=0;}
    if (k0 < noutedge || !nwp->directed_flag){
      for(e = EdgetreeMinimum(nwp->outedges, head);
      (fvalid==1 && ((tail1 = nwp->outedges[e].value) != 0));
      e = EdgetreeSuccessor(nwp->outedges, e)){
        if(alter==tail1){fvalid=0;}}
    }
    if (k0 >= noutedge || !nwp->directed_flag){
      for(e = EdgetreeMinimum(nwp->inedges, head);
      (fvalid==1 && ((tail1 = nwp->inedges[e].value) != 0));
      e = EdgetreeSuccessor(nwp->inedges, e)){
        if(alter==tail1){fvalid=0;}}
    }
    k1++;
  }

  if (k1 == 100){
    fvalid=0;
    continue;
  }
  
  if ( (!nwp->directed_flag && alter > head) ||
       (nwp->directed_flag && k0 < noutedge) )
    {
      MHp->togglehead[1] = head;
      MHp->toggletail[1] = alter;
    }else{
      MHp->togglehead[1] = alter;
      MHp->toggletail[1] = head;
    }
  
  if (!nwp->directed_flag){
    /* Check undirected degrees */
    k0 =nwp->outdegree[head]  + nwp->indegree[head];
    j0h=nwp->outdegree[tail]  + nwp->indegree[tail];
    j1h=nwp->outdegree[alter] + nwp->indegree[alter];
    
    j0=j0h-1;
    j1=j1h+1;
    
    if( ( (j0==j1h) && (j1==j0h) ) ){
      fvalid = 1;
    }else{
      fvalid = 0;
    }
  }else{
    /* Check directed degrees */
   if(k0 < noutedge){
     /* Check indegrees */
     j0h=nwp->indegree[tail];
     j1h=nwp->indegree[alter];
   }else{
     /* Check outdegrees */
     j0h=nwp->outdegree[tail];
     j1h=nwp->outdegree[alter];
   }
   j0=j0h-1;
   j1=j1h+1;
   
   if( ( (j0==j1h) && (j1==j0h) ) ){
     fvalid = 1;
   }else{
     fvalid = 0;
   }
  }
  
  }

  if (trynode==500){
    MHp->togglehead[1] = MHp->togglehead[0];
    MHp->toggletail[1] = MHp->toggletail[0];
  }
}

/*********************
 void MH_CondOutDegreeDist
*********************/
void MH_CondOutDegreeDist (MHproposal *MHp, DegreeBound *bd, Network *nwp) {  
  int noutedge=0, k, fvalid=0;
  int k0, k1;
  int trynode;
  Vertex e, alter, head=0, tail, tail1;
  MHp->ratio=1.0;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=2;    
    return;
  }

  fvalid = 0;
  trynode = 0;
  while(fvalid==0 && trynode < 1500){

  trynode++;

  while(noutedge==0){
    /* select a node at random */
    head = 1 + unif_rand() * nwp->nnodes;
    noutedge = nwp->outdegree[head];
  }
  
  k0 = (int)(unif_rand() * noutedge); 
  k=0;
  for(e = EdgetreeMinimum(nwp->outedges, head);
      ((tail = nwp->outedges[e].value) != 0 && k<k0);
      e = EdgetreeSuccessor(nwp->outedges, e)){++k;}
  MHp->togglehead[0] = head;
  MHp->toggletail[0] = tail;
  
  k1=0;
  fvalid=0;
  while(fvalid==0 && k1 < 100){
    while((alter = 1 + unif_rand() * nwp->nnodes) == head);
    fvalid=1;
    if(alter == tail){fvalid=0;}
    for(e = EdgetreeMinimum(nwp->outedges, head);
	(fvalid==1 && ((tail1 = nwp->outedges[e].value) != 0));
	e = EdgetreeSuccessor(nwp->outedges, e)){
      if(alter==tail1){fvalid=0;}}
    k1++;
  }
  if (k1 == 100){
    fvalid=0;
    continue;
  }
  
  MHp->togglehead[1] = head;
  MHp->toggletail[1] = alter;
  }
  
  if(trynode==1500 || !CheckTogglesValid(MHp, bd, nwp)){
      MHp->togglehead[0] = 1;
      MHp->toggletail[0] = 2;
      MHp->togglehead[1] = 1;
      MHp->toggletail[1] = 2;
  }
  

}

/*********************
 void MH_CondInDegreeDist
*********************/
void MH_CondInDegreeDist (MHproposal *MHp, DegreeBound *bd, Network *nwp) {  
  int ninedge=0, k, fvalid=0;
  int k0, k1;
  int trynode;
  Vertex e, alter, head=0, tail, tail1;
  MHp->ratio=1.0;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=2;    
    return;
  }

  fvalid = 0;
  trynode = 0;
  while(fvalid==0 && trynode < 1500){

  trynode++;

  while(ninedge==0){
    /* select a node at random */
    head = 1 + unif_rand() * nwp->nnodes;
    ninedge = nwp->indegree[head];
  }
  
  k0 = (int)(unif_rand() * ninedge); 
  k=0;
  for(e = EdgetreeMinimum(nwp->inedges, head);
      ((tail = nwp->inedges[e].value) != 0 && k<k0);
      e = EdgetreeSuccessor(nwp->inedges, e)){++k;}
  MHp->togglehead[0] = tail;
  MHp->toggletail[0] = head;
  
  k1=0;
  fvalid=0;
  while(fvalid==0 && k1 < 100){
    while((alter = 1 + unif_rand() * nwp->nnodes) == head);
    fvalid=1;
    if(alter == tail){fvalid=0;}
    for(e = EdgetreeMinimum(nwp->inedges, head);
	(fvalid==1 && ((tail1 = nwp->inedges[e].value) != 0));
	e = EdgetreeSuccessor(nwp->inedges, e)){
      if(alter==tail1){fvalid=0;}}
    k1++;
  }
  if (k1 == 100){
    fvalid=0;
    continue;
  }
  
  MHp->togglehead[1] = alter;
  MHp->toggletail[1] = head;
  
  }
  
  if(trynode==1500){
      MHp->togglehead[0] = 1;
      MHp->toggletail[0] = 2;
      MHp->togglehead[1] = 1;
      MHp->toggletail[1] = 2;
  }
}

/*********************
 void MH_CondDegree
*********************/
void MH_CondDegree (MHproposal *MHp, DegreeBound *bd, Network *nwp)  {  
  MHp->ratio=1.0;
  
  if(MHp->ntoggles == 0) { /* Initialize CondDeg by */
	                   /* Choosing Hexad or Tetrad */
    if( unif_rand() > 0.9 ){
      MHp->ntoggles=6;
    }else{
      MHp->ntoggles=4;
    }
    return;
  }

  if(MHp->ntoggles == 6) { /* Call Hexad */
    MH_CondDegHexadToggles (MHp, bd, nwp);
  }else{ /* call Tetrad */
    MH_CondDegTetradToggles (MHp, bd, nwp);
  }
}

/*********************
 void MH_CondDegHexadToggles
*********************/
void MH_CondDegHexadToggles (MHproposal *MHp, DegreeBound *bd, Network *nwp)  {  
  int x1, x2, x3, x4, x5, x6;
  int fvalid, trynode;
  Vertex head1, head2, head3, tail1, tail2, tail3;
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
  
  head1 = 1 + unif_rand() * nwp->nnodes;
  while((tail2 = 1 + unif_rand() * nwp->nnodes) == head1);
  tail3 = 1 + unif_rand() * nwp->nnodes;
  while(tail3 == tail2 || tail3 == head1){
    tail3 = 1 + unif_rand() * nwp->nnodes;
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
    tail1 = 1 + unif_rand() * nwp->nnodes;
    while(tail1 == tail2 || tail1 == tail3 || tail1 == head1){
      tail1 = 1 + unif_rand() * nwp->nnodes;
    }
    head2 = 1 + unif_rand() * nwp->nnodes;
    while(head2 == tail2 || head2 == tail3 || head2 == head1 ||
	  head2 == tail1 ){
      head2 = 1 + unif_rand() * nwp->nnodes;
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
	head3 = 1 + unif_rand() * nwp->nnodes;
	while(head3 == tail2 || head3 == tail3 || head3 == head1 ||
	      head3 == tail1 || head3 == head2 ){
	  head3 = 1 + unif_rand() * nwp->nnodes;
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

/*********************
 void MH_CondDegTetradToggles
*********************/
void MH_CondDegTetradToggles (MHproposal *MHp, DegreeBound *bd, Network *nwp) {  
  int x1, x2, x3, x4;
  int fvalid, trynode;
  Vertex head1, head2, tail1, tail2;
  MHp->ratio=1.0;
  
  fvalid = 0;
  trynode = 0;
  while(fvalid==0 && trynode < 5000){

  trynode++;
  /* select a node at random */
  head1 = 1 + unif_rand() * nwp->nnodes;
  while((tail1 = 1 + unif_rand() * nwp->nnodes) == head1);
  tail2 = 1 + unif_rand() * nwp->nnodes;
  while(tail2 == tail1 || tail2 == head1){
    tail2 = 1 + unif_rand() * nwp->nnodes;
  }
  if ((!nwp->directed_flag) && head1 > tail1){
    x1 = EdgetreeSearch(tail1, head1, nwp->outedges) > 0;
  }else{
    x1 = EdgetreeSearch(head1, tail1, nwp->outedges) > 0;
  }
  if ((!nwp->directed_flag) && head1 > tail2){
    x2 = EdgetreeSearch(tail2, head1, nwp->outedges) > 0;
  }else{
    x2 = EdgetreeSearch(head1, tail2, nwp->outedges) > 0;
  }
  if (x1 != x2){
    head2 = 1 + unif_rand() * nwp->nnodes;
    while(head2 == tail2 || head2 == tail1 || head2 == head1 ){
      head2 = 1 + unif_rand() * nwp->nnodes;
    }
    if ((!nwp->directed_flag) && head2 > tail1){
      x3 = EdgetreeSearch(tail1, head2, nwp->outedges) > 0;
    }else{
      x3 = EdgetreeSearch(head2, tail1, nwp->outedges) > 0;
    }
    if (x2 == x3){
      if ((!nwp->directed_flag) && head2 > tail2){
	x4 = EdgetreeSearch(tail2, head2, nwp->outedges) > 0;
      }else{
	x4 = EdgetreeSearch(head2, tail2, nwp->outedges) > 0;
      }
      if (x4 == x1){
	if ( (!nwp->directed_flag) ){
	  if ( head1 > tail1 ){
	    MHp->togglehead[0] = tail1;
	    MHp->toggletail[0] = head1;
	  }else{
	    MHp->togglehead[0] = head1;
	    MHp->toggletail[0] = tail1;
	  }
	  if ( head1 > tail2 ){
	    MHp->togglehead[1] = tail2;
	    MHp->toggletail[1] = head1;
	  }else{
	    MHp->togglehead[1] = head1;
	    MHp->toggletail[1] = tail2;
	  }
	  if ( head2 > tail1 ){
	    MHp->togglehead[2] = tail1;
	    MHp->toggletail[2] = head2;
	  }else{
	    MHp->togglehead[2] = head2;
	    MHp->toggletail[2] = tail1;
	  }
	  if ( head2 > tail2 ){
	    MHp->togglehead[3] = tail2;
	    MHp->toggletail[3] = head2;
	  }else{
	    MHp->togglehead[3] = head2;
	    MHp->toggletail[3] = tail2;
	  }
	}else{
	  MHp->togglehead[0] = head1;
	  MHp->toggletail[0] = tail1;
	  MHp->togglehead[1] = head1;
	  MHp->toggletail[1] = tail2;
	  MHp->togglehead[2] = head2;
	  MHp->toggletail[2] = tail1;
	  MHp->togglehead[3] = head2;
	  MHp->toggletail[3] = tail2;
	}
	fvalid=1;
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
  }
}
/*********************
 void MH_TwoRandomToggles
*********************/
void MH_TwoRandomToggles (MHproposal *MHp, DegreeBound *bd, Network *nwp) {  
  Vertex head, tail;
  int i;
  
  if(MHp->ntoggles == 0) { /* Initialize OneRandomToggle */
    MHp->ntoggles=2;
    return;
  }
  MHp->ratio = 1.0;

  for (i = 0; i < 2; i++){
   head = 1 + unif_rand() * nwp->nnodes;
   while ((tail = 1 + unif_rand() * nwp->nnodes) == head);
   if (!nwp->directed_flag && head > tail) {
     MHp->togglehead[i] = tail;
     MHp->toggletail[i] = head;
   }else{
     MHp->togglehead[i] = head;
     MHp->toggletail[i] = tail;
   }
  }
}

/*********************
 void MH_RandomNode
*********************/
void MH_randomnode (MHproposal *MHp, DegreeBound *bd, Network *nwp) {
  
  Vertex root, alter;
  int j;
  
  if(MHp->ntoggles == 0) { /* Initialize OneRandomToggle */
    MHp->ntoggles= nwp->nnodes - 1;
    return;
  }
  MHp->ratio = 1.0;

  root = 1 + unif_rand() * nwp->nnodes;
  
  j = 0;
  for (alter = 1; alter <= nwp->nnodes; alter++)
    {
      /* there is never an edge (root, root) */
      if (alter != root) {
       if (!nwp->directed_flag && root > alter) {
        MHp->togglehead[j] = alter;
        MHp->toggletail[j] = root;
       }else{
        MHp->togglehead[j] = root;
        MHp->toggletail[j] = alter;
       }
       j++;
      }
    }
}

void MH_randomtoggleNonObserved (MHproposal *MHp, DegreeBound *bd, Network *nwp)  {  
  Edge rane, nmissing = nwp[1].nedges;
  
  if(MHp->ntoggles == 0) { /* Initialize randomtoggle */
    MHp->ntoggles=1;
    return;
  }
  MHp->ratio = 1.0;

  if(nmissing==0){
    *MHp->togglehead = MH_FAILED;
    *MHp->toggletail = MH_IMPOSSIBLE;
  }

  rane = 1 + unif_rand() * nmissing;
  FindithEdge(MHp->togglehead, MHp->toggletail, rane, &nwp[1]);
}

/* The ones below have not been tested */

/*********************
 void MH_ConstrainedCondOutDegDist
*********************/
void MH_ConstrainedCondOutDegDist (MHproposal *MHp, DegreeBound *bd, Network *nwp){  
  int noutedge=0, k, fvalid=0;
  int k0, k1;
  Vertex e, alter, head, tail, tail1;
  MHp->ratio=1.0;
  
  while(noutedge==0){
    /* select a node at random */
    head = 1 + unif_rand() * nwp->nnodes;
    noutedge = nwp->outdegree[head];
  }
  
  k0 = (int)(unif_rand() * noutedge); 
  k=0;
  for(e = EdgetreeMinimum(nwp->outedges, head);
      ((tail = nwp->outedges[e].value) != 0 && k<k0);
      e = EdgetreeSuccessor(nwp->outedges, e)){++k;}
  MHp->togglehead[0] = head;
  MHp->toggletail[0] = tail;
  
  k1=0;
  fvalid=0;
  while(fvalid==0 && k1 < 100){
    while((alter = 1 + unif_rand() * nwp->nnodes) == head);
    fvalid=1;
    if(alter == tail){fvalid=0;}
    for(e = EdgetreeMinimum(nwp->outedges, head);
	(fvalid==1 && ((tail1 = nwp->outedges[e].value) != 0));
	e = EdgetreeSuccessor(nwp->outedges, e)){
      if(alter==tail1){fvalid=0;}}
    k1++;
  }
  if (k1 == 100){
    MHp->togglehead[0] = MHp->toggletail[0] = 0;
    MHp->togglehead[1] = MHp->toggletail[1] = 0;
  }
  
  MHp->togglehead[1] = head;
  MHp->toggletail[1] = alter;
  
  if (!fvalid){
    MHp->togglehead[0] = MHp->toggletail[0] = 0;
    MHp->togglehead[1] = MHp->toggletail[1] = 0;
  }
  
  for(k=0; k < 2; k++){
    if (DesignMissing(MHp->togglehead[k], MHp->toggletail[k], &nwp[1])==0){
      MHp->togglehead[0] = MHp->toggletail[0] = 0;
      MHp->togglehead[1] = MHp->toggletail[1] = 0;
    }
  }
}


void MH_NodePairedTiesToggles (MHproposal *MHp, DegreeBound *bd, Network *nwp) {  
  /* chooses a node and toggles all ties and
	 and toggles an equal number of matching nonties
	 for that node */
  int nedge=0,j,k;
  int fvalid = 1;
  Vertex e, head, prop;
  MHp->ratio=1.0;
  
  /* double to integer coercion */
  head = 1 + unif_rand() * nwp->nnodes; 
  
  for(e = EdgetreeMinimum(nwp->outedges, head);
      (prop = nwp->outedges[e].value) != 0; /* loop if */
      e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of head */
    {
      MHp->togglehead[nedge] = head;
      MHp->toggletail[nedge] = prop;
      ++nedge;
    }
  for(e = EdgetreeMinimum(nwp->inedges, head);
      (prop = nwp->inedges[e].value) != 0; /* loop if */
      e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of head */
    {
      MHp->toggletail[nedge] = head;
      MHp->togglehead[nedge] = prop;
      ++nedge;
    }
  
  if(nedge > nwp->nnodes-nedge){
    MHp->togglehead[0] = MHp->toggletail[0] = 0;
    MHp->togglehead[1] = MHp->toggletail[1] = 0;
  }  
  j = 0;
  while (j <=nedge)
    {
      prop = 1 + unif_rand() * nwp->nnodes; 
      k=0;
      fvalid=1;
      while(fvalid==1 && k<nedge+j){
	if(EdgetreeSearch( MIN(prop,MHp->togglehead[k]),
			   MAX(prop,MHp->togglehead[k]), nwp->outedges) +
	   EdgetreeSearch( MIN(prop,MHp->toggletail[k]),
			   MAX(prop,MHp->toggletail[k]), nwp->outedges)==0
	   ){++k;
	}else{
	  fvalid=0;
	}
      }
      if(prop>head){
	MHp->togglehead[j+nedge] = head;
	MHp->toggletail[j+nedge] = prop;
      }else{
	MHp->togglehead[j+nedge] = prop;
	MHp->toggletail[j+nedge] = head;
      }
      ++j;
    }
  
  j = 2*nedge;
  if (!CheckTogglesValid(MHp, bd, nwp))
    {
      *MHp->togglehead = *MHp->toggletail = 0;
    }
}

/*********************
 void MH_OneRandomTnTNode
*********************/
void MH_OneRandomTnTNode (MHproposal *MHp, DegreeBound *bd, Network *nwp) {  
  Vertex head=0, tail, e, tail1;
  int noutedge=0, ninedge=0, k0=0, ndyad, fvalid=0, k;
  
  if ( nwp->directed_flag )
    {
      ndyad = (nwp->nnodes - 1) * nwp->nnodes;
    }else{
      ndyad = (nwp->nnodes - 1) * nwp->nnodes / 2;
    }
  
  fvalid=0;
  while(fvalid==0){
    
    if ( unif_rand() < 0.5 && nwp->nedges > 0) 
      {
	
	/* select a tie */
	ninedge=0;
	noutedge=0;
	while(noutedge+ninedge==0){
	  /* select a node at random */
	  head = 1 + unif_rand() * nwp->nnodes;
	  ninedge = nwp->indegree[head];
	  noutedge = nwp->outdegree[head];
	}
	
	k0 = (int)(unif_rand() * (noutedge+ninedge)); 
	if (k0 < noutedge){
	  k=0;
	  for(e = EdgetreeMinimum(nwp->outedges, head);
	      ((tail = nwp->outedges[e].value) != 0 && k<k0);
	      e = EdgetreeSuccessor(nwp->outedges, e)){++k;}
	}else{
	  k=0;
	  for(e = EdgetreeMinimum(nwp->inedges, head);
	      ((tail = nwp->inedges[e].value) != 0 && k<(k0-noutedge));
	      e = EdgetreeSuccessor(nwp->inedges, e)){++k;}
	}
	if ( (!nwp->directed_flag && head > tail) ||
	     (nwp->directed_flag && k0 >= noutedge) )
	  {
	    MHp->togglehead[0] = tail;
	    MHp->toggletail[0] = head;
	  }else{
	    MHp->togglehead[0] = head;
	    MHp->toggletail[0] = tail;
	  }
	
	MHp->ratio = ((noutedge+ninedge)*1.0)/(nwp->nnodes-1-noutedge-ninedge-1);
	fvalid =1;
      }else{
	/* Choose random non-tie */

	/* select a node at random */
	ninedge=nwp->nnodes-1;
	noutedge=0;
	while(noutedge+ninedge>=(nwp->nnodes-1)){
	  ninedge=0;
	  /* select a node at random */
	  head = 1 + unif_rand() * nwp->nnodes;
	  ninedge = nwp->indegree[head];
	  noutedge = nwp->outdegree[head];
	}
	
	fvalid=0;
	while(fvalid==0){
	  while ((tail = 1 + unif_rand() * nwp->nnodes) == head);
	  fvalid=1;
	  for(e = EdgetreeMinimum(nwp->outedges, head);
	      (fvalid==1 && ((tail1 = nwp->outedges[e].value) != 0));
	      e = EdgetreeSuccessor(nwp->outedges, e)){
	    if(tail==tail1){fvalid=0;}}
	  if (!(nwp->directed_flag)){
	    for(e = EdgetreeMinimum(nwp->inedges, head);
		(fvalid==1 && ((tail1 = nwp->inedges[e].value) != 0));
		e = EdgetreeSuccessor(nwp->inedges, e)){
	      if(tail==tail1){fvalid=0;}}
	  }
	}
	
	if ( (!nwp->directed_flag && head > tail) ||
	     (nwp->directed_flag && k0 >= noutedge) )
	  {
	    MHp->togglehead[0] = tail;
	    MHp->toggletail[0] = head;
	  }else{
	    MHp->togglehead[0] = head;
	    MHp->toggletail[0] = tail;
	  }
	
        if ( nwp->directed_flag )
	  {
	    MHp->ratio = (nwp->nnodes-1-noutedge-ninedge)/(noutedge+ninedge+1.0);
	  }else{
	    MHp->ratio = (nwp->nnodes-1-noutedge-ninedge)/(noutedge+ninedge+1.0);
	  }
      }
  }
}

/*********************
 void MH_ReallocateWithReplacement
*********************/
void MH_ReallocateWithReplacement (MHproposal *MHp, DegreeBound *bd, Network *nwp) {  
  int i;
  Vertex root;
  Vertex* edges;
  int edgecount = 0;
  
  MHp->ratio=1.0;
  /* select a node at random */
  root = 1 + unif_rand() * nwp->nnodes;

  edges = (Vertex *) malloc(sizeof(Vertex) * (nwp->nnodes+1));
  for (i = 0; i <= nwp->nnodes; i++)
    edges[i] = NO_EDGE;
  
  /* count current edges and mark them in an array */
  for (i = 1; i <= nwp->nnodes; i++)
    {
      if (root == i) continue;
      if (EdgetreeSearch(root, i, nwp->outedges) > 0)
	{
	  edges[i] = OLD_EDGE;
	  edgecount++;
	}
      if (!nwp->directed_flag && (root > i) &&
	  (EdgetreeSearch(i, root, nwp->outedges) > 0))
	{
	  edges[i] = OLD_EDGE;
	  edgecount++;
	}
    }
  
  /* select edgecount edges to create */
  for (i = 0; i < edgecount; i++)
    {
      Vertex newtail;
      /* get a new edge, neither the root nor something already chosen */
      while ((newtail = 1 + unif_rand() * nwp->nnodes) == root ||
	     (edges[newtail] & NEW_EDGE))
	;
      
      /* if this edge already exists - (OLD_EDGE | NEW_EDGE) == CAN_IGNORE */
      edges[newtail] = edges[newtail] | NEW_EDGE;
    }
  
  /* index into MHp->togglehead/MHp->toggletail is  */
  edgecount = 0;
  
  /* add to toggle list:  anything that is non zero in edges array
     should be toggled, whether on or off. */
  for (i = 0; i <= nwp->nnodes; i++)
    {
      if (edges[i] == NO_EDGE || edges[i] == CAN_IGNORE) continue;
      
      /* double to integer coercion */
      MHp->togglehead[edgecount] = root;
      MHp->toggletail[edgecount] = i;
      
      if (!nwp->directed_flag && (MHp->togglehead[edgecount] > MHp->toggletail[edgecount]))
	{
	  Vertex temp;
	  temp = MHp->togglehead[edgecount];
	  MHp->togglehead[edgecount] = MHp->toggletail[edgecount];
	  MHp->toggletail[edgecount] = temp;
	}
      edgecount++;
    }
  free(edges);
}

/*********************
 void MH_AllTogglesForOneNode
*********************/
void MH_AllTogglesForOneNode (MHproposal *MHp, DegreeBound *bd, Network *nwp) {
  
  int i;
  int j;
  int root;
  
  MHp->ratio=1.0;
  root = 1 + unif_rand() * nwp->nnodes;
  
  j = 0;
  for (i = 1; i <= nwp->nnodes; i++)
    {
      /* probability here only do this with .8? */
      
      /* there is never an edge (root, root) */
      if (i == root)
	continue;
      
      /* double to integer coercion */
      MHp->togglehead[j] = root;
      MHp->toggletail[j] = i;
      
      if (!nwp->directed_flag && (MHp->togglehead[j] > MHp->toggletail[j]))
	{
	  Vertex temp;
	  temp = MHp->togglehead[j];
	  MHp->togglehead[j] = MHp->toggletail[j];
	  MHp->toggletail[j] = temp;
	}
      j++;
    }
}


/*********************
 void MH_SwitchLabelTwoNodesToggles
*********************/
void MH_SwitchLabelTwoNodesToggles (MHproposal *MHp, DegreeBound *bd, Network *nwp) {  
  int nedge1=0, nedge2=0, k, ntoggles;
  Vertex *edges1, *edges2;
  Vertex e, head2, tail2, head1, tail1;
  MHp->ratio=1.0;
  
  /* select a node at random */
  edges1 = (Vertex *) malloc(sizeof(Vertex) * (nwp->nnodes+1));
  edges2 = (Vertex *) malloc(sizeof(Vertex) * (nwp->nnodes+1));
  
  while(nedge1==0){
    head1 = 1 + unif_rand() * nwp->nnodes;
    
    for(e = EdgetreeMinimum(nwp->outedges, head1);
	(tail1 = nwp->outedges[e].value) != 0; /* loop if */
	e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of head1 */
      {
        edges1[nedge1] = tail1;
	++nedge1;
      }
    for(e = EdgetreeMinimum(nwp->inedges, head1);
	(tail1 = nwp->inedges[e].value) != 0; /* loop if */
	e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of head1 */
      {
        edges1[nedge1] = tail1;
	++nedge1;
      }
  }
  
  while((head2 = 1 + unif_rand() * nwp->nnodes) == head1);
  
  for(e = EdgetreeMinimum(nwp->outedges, head2);
      (tail2 = nwp->outedges[e].value) != 0; /* loop if */
      e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of head2 */
    {
      edges2[nedge2] = tail2;
      ++nedge2;
    }
  for(e = EdgetreeMinimum(nwp->inedges, head2);
      (tail2 = nwp->inedges[e].value) != 0; /* loop if */
      e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of head2 */
    {
      edges2[nedge2] = tail2;
      ++nedge2;
    }
  
  ntoggles = 0;
  for(k=0; k < nedge1; k++){
    if (head1 > edges1[k])
      {
	MHp->togglehead[ntoggles] = edges1[k];
	MHp->toggletail[ntoggles] = head1;
      }
    if (head1 < edges1[k]){
      MHp->togglehead[ntoggles] = head1;
      MHp->toggletail[ntoggles] = edges1[k];
    }
    if(head1 != edges1[k]) ntoggles++;
  }
  
  for(k=0; k < nedge2; k++){
    if (head1 > edges2[k])
      {
	MHp->togglehead[ntoggles] = edges2[k];
	MHp->toggletail[ntoggles] = head1;
      }
    if (head1 < edges2[k]){
      MHp->togglehead[ntoggles] = head1;
      MHp->toggletail[ntoggles] = edges2[k];
    }
    if(head1 != edges2[k]) ntoggles++;
  }
  
  for(k=0; k < nedge2; k++){
    if (head2 > edges2[k])
      {
	MHp->togglehead[ntoggles] = edges2[k];
	MHp->toggletail[ntoggles] = head2;
      }
    if (head2 < edges2[k]){
      MHp->togglehead[ntoggles] = head2;
      MHp->toggletail[ntoggles] = edges2[k];
    }
    if(head2 != edges2[k]) ntoggles++;
  }
  
  for(k=0; k < nedge1; k++){
    if (head2 > edges1[k])
      {
	MHp->togglehead[ntoggles] = edges1[k];
	MHp->toggletail[ntoggles] = head2;
      }
    if (head2 < edges1[k]){
      MHp->togglehead[ntoggles] = head2;
      MHp->toggletail[ntoggles] = edges1[k];
    }
    if(head2 != edges1[k]) ntoggles++;
  }
  free(edges1);
  free(edges2);
}


/*********************
 void MH_ConstrainedCondDegDist
*********************/
void MH_ConstrainedCondDegDist (MHproposal *MHp, DegreeBound *bd, Network *nwp)  {  
  int noutedge=0, ninedge=0, k, fvalid=0;
  int k0, j0, j1, k1;
  int j0h, j1h;
  Vertex *outedges, *inedges;
  Vertex e, alter, head=0, tail;
  MHp->ratio=1.0;
  
  /* select a node at random */
  outedges = (Vertex *) malloc(sizeof(Vertex) * (nwp->nnodes+1));
  inedges = (Vertex *) malloc(sizeof(Vertex) * (nwp->nnodes+1));
  
  while(noutedge==0 && ninedge==0){
    head = 1 + unif_rand() * nwp->nnodes;
    
    for(e = EdgetreeMinimum(nwp->outedges, head);
	(tail = nwp->outedges[e].value) != 0; /* loop if */
	e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of head */
      {
        outedges[noutedge] = tail;
	++noutedge;
      }
    for(e = EdgetreeMinimum(nwp->inedges, head);
	(tail = nwp->inedges[e].value) != 0; /* loop if */
	e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of head */
      {
        inedges[ninedge] = tail;
	++ninedge;
      }
  }
  
  k0 = (int)(unif_rand() * (noutedge+ninedge)); 
  if (k0 < noutedge){
    tail = outedges[k0]; 
  }else{
    tail = inedges[k0-noutedge]; 
  }
  if ( (!nwp->directed_flag && head > tail) ||
       (  nwp->directed_flag  && k0 >= noutedge) )
    {
      MHp->togglehead[0] = tail;
      MHp->toggletail[0] = head;
    }else{
      MHp->togglehead[0] = head;
      MHp->toggletail[0] = tail;
    }
  
  if (DesignMissing(MHp->togglehead[0], MHp->toggletail[0], &nwp[1])==0){
    MHp->togglehead[0] = MHp->toggletail[0] = 0;
    MHp->togglehead[1] = MHp->toggletail[1] = 0;
  }
  
  fvalid=0;
  k1=0;
  while(fvalid==0 && k1 < 100){
    while((alter = 1 + unif_rand() * nwp->nnodes) == head);
    if(alter != tail){fvalid=1;}
    fvalid=1;
    if (k0 < noutedge || !(nwp->directed_flag)){
      k=0;
      while(fvalid==1 && noutedge > 0 && k <= noutedge-1){
	if(alter == outedges[k]){fvalid=0;}else{++k;}
      }
    }
    if (k0 >= noutedge || !(nwp->directed_flag)){
      k=0;
      while(fvalid==1 && ninedge > 0 && k <= ninedge-1){
	if(alter == inedges[k]){fvalid=0;}else{++k;}
      }
    }
    k1++;
  }
  
  if (k1 == 100){
    MHp->togglehead[0] = MHp->toggletail[0] = 0;
    MHp->togglehead[1] = MHp->toggletail[1] = 0;
  }
  
  if ( (!nwp->directed_flag && alter > head) ||
       (nwp->directed_flag && k0 < noutedge) )
    {
      MHp->togglehead[1] = head;
      MHp->toggletail[1] = alter;
    }else{
      MHp->togglehead[1] = alter;
      MHp->toggletail[1] = head;
    }
  
  if (DesignMissing(MHp->togglehead[1], MHp->toggletail[1], &nwp[1])==0){
    MHp->togglehead[0] = MHp->toggletail[0] = 0;
    MHp->togglehead[1] = MHp->toggletail[1] = 0;
  }
  
  free(outedges);
  free(inedges);
  
  /* Check undirected degrees */
  if (!nwp->directed_flag){
    k0=nwp->outdegree[head]+ nwp->indegree[head];
    j0h=nwp->outdegree[tail]+ nwp->indegree[tail];
    j1h=nwp->outdegree[alter]+ nwp->indegree[alter];
    
    j0=j0h-1;
    j1=j1h+1;
    
    if( ( (j0==j1h) && (j1==j0h) ) ){
      fvalid = 1;
    }else{
      fvalid = 0;
    }
  }else{
    if(k0 < noutedge){
      /* Check indegrees */
      j0h=nwp->indegree[tail];
      j1h=nwp->indegree[alter];
    }else{
      /* Check outdegrees */
      j0h=nwp->outdegree[tail];
      j1h=nwp->outdegree[alter];
    }
    j0=j0h-1;
    j1=j1h+1;
    
    if( ( (j0==j1h) && (j1==j0h) ) ){
      fvalid = 1;
    }else{
      fvalid = 0;
    }
  }
  
  if (!fvalid){
    MHp->togglehead[0] = MHp->toggletail[0] = 0;
    MHp->togglehead[1] = MHp->toggletail[1] = 0;
  }
}

void MH_ConstrainedNodePairedTiesToggles (MHproposal *MHp,
       	       DegreeBound *bd, Network *nwp) {  
  /* chooses a node and toggles all ties and
     and toggles an equal number of matching nonties
     for that node */
  int nedge=0,j,k;
  int fvalid = 1;
  Vertex e, head, prop;
  MHp->ratio=1.0;
  
  /* double to integer coercion */
  head = 1 + unif_rand() * nwp->nnodes; 
  
  for(e = EdgetreeMinimum(nwp->outedges, head);
      (prop = nwp->outedges[e].value) != 0; /* loop if */
      e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of head */
    {
      MHp->togglehead[nedge] = head;
      MHp->toggletail[nedge] = prop;
      ++nedge;
    }
  for(e = EdgetreeMinimum(nwp->inedges, head);
      (prop = nwp->inedges[e].value) != 0; /* loop if */
      e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of head */
    {
      MHp->toggletail[nedge] = head;
      MHp->togglehead[nedge] = prop;
      ++nedge;
    }
  
  if(nedge > nwp->nnodes-nedge){
    MHp->togglehead[0] = MHp->toggletail[0] = 0;
    MHp->togglehead[1] = MHp->toggletail[1] = 0;
  }  
  j = 0;
  while (j <=nedge)
    {
      prop = 1 + unif_rand() * nwp->nnodes; 
      k=0;
      fvalid=1;
      while(fvalid==1 && k<nedge+j){
	if(EdgetreeSearch( MIN(prop,MHp->togglehead[k]),
			   MAX(prop,MHp->togglehead[k]), nwp->outedges) +
	   EdgetreeSearch( MIN(prop,MHp->toggletail[k]),
			   MAX(prop,MHp->toggletail[k]), nwp->outedges)==0
	   ){++k;
	}else{
	  fvalid=0;}
      }
      if(prop>head){
	MHp->togglehead[j+nedge] = head;
	MHp->toggletail[j+nedge] = prop;
      }else{
	MHp->togglehead[j+nedge] = prop;
	MHp->toggletail[j+nedge] = head;
      }
      ++j;
    }
  
  j = 2*nedge;
  if (!CheckConstrainedTogglesValid(MHp, bd, nwp))
    {
      *MHp->togglehead = *MHp->toggletail = 0;
    }
}

/*********************
 void MH_ConstrainedReallocateWithReplacement
*********************/
void MH_ConstrainedReallocateWithReplacement (MHproposal *MHp,
       	       DegreeBound *bd, Network *nwp) {  
  int i;
  Vertex root;
  Vertex* edges;
  int edgecount = 0;
  
  MHp->ratio=1.0;
  /* select a node at random */
  root = 1 + unif_rand() * nwp->nnodes;

  edges = (Vertex *) malloc(sizeof(Vertex) * (nwp->nnodes+1));
  for (i = 0; i <= nwp->nnodes; i++)
    edges[i] = NO_EDGE;
  
  /* count current edges and mark them in an array */
  for (i = 1; i <= nwp->nnodes; i++)
    {
      if (root == i) continue;
      if (EdgetreeSearch(root, i, nwp->outedges) > 0)
	{
	  edges[i] = OLD_EDGE;
	  edgecount++;
	}
      if (!nwp->directed_flag && (root > i) &&
	  (EdgetreeSearch(i, root, nwp->outedges) > 0))
	{
	  edges[i] = OLD_EDGE;
	  edgecount++;
	}
    }
  
  /* select edgecount edges to create */
  for (i = 0; i < edgecount; i++)
    {
      Vertex newtail;
      
      /* get a new edge, neither the root nor something already chosen */
      while ((newtail = 1 + unif_rand() * nwp->nnodes) == root ||
	     (edges[newtail] & NEW_EDGE))
	;
      
      /* if this edge already exists - (OLD_EDGE | NEW_EDGE) == CAN_IGNORE */
      edges[newtail] = edges[newtail] | NEW_EDGE;
    }
  
  /* index into MHp->togglehead/MHp->toggletail is  */
  edgecount = 0;
  
  /* add to toggle list:  anything that is non zero in edges array
     should be toggled, whether on or off. */
  for (i = 0; i <= nwp->nnodes; i++)
    {
      if (edges[i] == NO_EDGE || edges[i] == CAN_IGNORE) continue;
      
      /* double to integer coercion */
      MHp->togglehead[edgecount] = root;
      MHp->toggletail[edgecount] = i;
      
      if (!nwp->directed_flag && (MHp->togglehead[edgecount] > MHp->toggletail[edgecount]))
	{
	  Vertex temp;
	  temp = MHp->togglehead[edgecount];
	  MHp->togglehead[edgecount] = MHp->toggletail[edgecount];
	  MHp->toggletail[edgecount] = temp;
	}
      edgecount++;
    }
  free(edges);
}

/*********************
 void MH_ConstrainedAllTogglesForOneNode
*********************/
void MH_ConstrainedAllTogglesForOneNode (MHproposal *MHp,
					 DegreeBound *bd, Network *nwp) {
  int i;
  int j;
  int root;
  
  MHp->ratio=1.0;
  root = 1 + unif_rand() * nwp->nnodes;
  
  j = 0;
  for (i = 1; i <= nwp->nnodes; i++)
    {
      /* probability here only do this with .8? */
      
      /* there is never an edge (root, root) */
      if (i == root)
	continue;
      
      /* double to integer coercion */
      MHp->togglehead[j] = root;
      MHp->toggletail[j] = i;
      
      if (!nwp->directed_flag && (MHp->togglehead[j] > MHp->toggletail[j]))
	{
	  Vertex temp;
	  temp = MHp->togglehead[j];
	  MHp->togglehead[j] = MHp->toggletail[j];
	  MHp->toggletail[j] = temp;
	}
      j++;
    }
}

/*********************
 void MH_ConstrainedTwoRandomToggles
*********************/
void MH_ConstrainedTwoRandomToggles (MHproposal *MHp,
				     DegreeBound *bd, Network *nwp) {  
  int i;
  MHp->ratio=1.0;
  
  for (i = 0; i < 2; i++)
    {
      /* double to integer coercion */
      MHp->togglehead[i] = 1 + unif_rand() * nwp->nnodes; 
      while ((MHp->toggletail[i] = 1 + unif_rand() * nwp->nnodes) == MHp->togglehead[i]);
      
      while(DesignMissing(MHp->togglehead[i], MHp->toggletail[i], &nwp[1])==0){
	MHp->togglehead[i] = 1 + unif_rand() * nwp->nnodes; 
	while ((MHp->toggletail[i] = 1 + unif_rand() * nwp->nnodes) == MHp->togglehead[i]);
      }
      if (!nwp->directed_flag && MHp->togglehead[i] > MHp->toggletail[i]) 
	{
	  Vertex temp;
	  temp = MHp->togglehead[i];
	  MHp->togglehead[i] = MHp->toggletail[i];
	  MHp->toggletail[i] = temp;
	}
    }
  
  if (!CheckConstrainedTogglesValid(MHp, bd, nwp))
    {
      MHp->togglehead[0] = MHp->toggletail[0] = 0;
      MHp->togglehead[1] = MHp->toggletail[1] = 0;
    }  
}

/*********************
 void MH_ConstrainedCondDeg
*********************/
void MH_ConstrainedCondDeg (MHproposal *MHp,
					 DegreeBound *bd, Network *nwp) {  
  /* WARNING: THIS NEEDS TO BE FIXED */
  int nedge1=0, nedge2=0, k, toomany, fvalid=0;
  Vertex *edges1, *edges2;
  Vertex e, head2=0, tail2, head1, tail1;
  MHp->ratio=1.0;
  
  /* select a node at random */
  edges1 = (Vertex *) malloc(sizeof(Vertex) * (nwp->nnodes+1));
  edges2 = (Vertex *) malloc(sizeof(Vertex) * (nwp->nnodes+1));
  
  while(nedge1==0){
    head1 = 1 + unif_rand() * nwp->nnodes;
    
    for(e = EdgetreeMinimum(nwp->outedges, head1);
	(tail1 = nwp->outedges[e].value) != 0; /* loop if */
	e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of head1 */
      {
        edges1[nedge1] = tail1;
	++nedge1;
      }
    for(e = EdgetreeMinimum(nwp->inedges, head1);
	(tail1 = nwp->inedges[e].value) != 0; /* loop if */
	e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of head1 */
      {
        edges1[nedge1] = tail1;
	++nedge1;
      }
  }
  
  tail1 = edges1[(int)(unif_rand() * nedge1)]; 
  if (head1 > tail1)
    {
      MHp->togglehead[0] = tail1;
      MHp->toggletail[0] = head1;
    }else{
      MHp->togglehead[0] = head1;
      MHp->toggletail[0] = tail1;
    }
   
  toomany = 0;
  while(nedge2==0 && toomany < 100){
    fvalid=0;
    while(fvalid==0){
      while((head2 = 1 + unif_rand() * nwp->nnodes) == head1);
      k=0;
      fvalid=1;
      while(fvalid==1 && k < nedge1){
	if(head2 == edges1[k]){fvalid=0;}else{++k;}
      }
    }

    for(e = EdgetreeMinimum(nwp->outedges, head2);
	(tail2 = nwp->outedges[e].value) != 0; /* loop if */
	e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of head2 */
      {
        edges2[nedge2] = tail2;
	++nedge2;
      }
    for(e = EdgetreeMinimum(nwp->inedges, head2);
	(tail2 = nwp->inedges[e].value) != 0; /* loop if */
	e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of head2 */
      {
        edges2[nedge2] = tail2;
	++nedge2;
      }
    ++toomany;
  }
  if (toomany==100){
    MHp->togglehead[0] = MHp->toggletail[0] = 0;
    MHp->togglehead[1] = MHp->toggletail[1] = 0;
  }
  toomany=0;
  fvalid=0;
  while(fvalid==0 && toomany < 10){
    while((tail2 = edges2[(int)(unif_rand() * nedge2)]) == head1);
    k=0;
    fvalid=1;
    while(fvalid==1 && k < nedge1){
      if(tail2 == edges1[k]){fvalid=0;}else{++k;}
    }
    ++toomany;
  }
  if (!fvalid || toomany==10){
    MHp->togglehead[0] = MHp->toggletail[0] = 0;
    MHp->togglehead[1] = MHp->toggletail[1] = 0;
    free(edges1);
    free(edges2);
      }
  if (head2 > tail2)
    {
      MHp->togglehead[1] = tail2;
      MHp->toggletail[1] = head2;
    }else{
      MHp->togglehead[1] = head2;
      MHp->toggletail[1] = tail2;
    }
  free(edges1);
  free(edges2);
}

/*********************
 void MH_ConstrainedSwitchLabelTwoNodesToggles
*********************/
void MH_ConstrainedSwitchLabelTwoNodesToggles (MHproposal *MHp,
       	       DegreeBound *bd, Network *nwp)  {  
  int nedge1=0, nedge2=0, k, ntoggles;
  Vertex *edges1, *edges2;
  Vertex e, head2, tail2, head1, tail1;
  MHp->ratio=1.0;
  
  /* select a node at random */

  edges1 = (Vertex *) malloc(sizeof(Vertex) * (nwp->nnodes+1));
  edges2 = (Vertex *) malloc(sizeof(Vertex) * (nwp->nnodes+1));

  while(nedge1==0){
    head1 = 1 + unif_rand() * nwp->nnodes;
    
    for(e = EdgetreeMinimum(nwp->outedges, head1);
	(tail1 = nwp->outedges[e].value) != 0; /* loop if */
	e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of head1 */
      {
        edges1[nedge1] = tail1;
	++nedge1;
      }
    for(e = EdgetreeMinimum(nwp->inedges, head1);
	(tail1 = nwp->inedges[e].value) != 0; /* loop if */
	e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of head1 */
      {
        edges1[nedge1] = tail1;
	++nedge1;
      }
  }
  
  while((head2 = 1 + unif_rand() * nwp->nnodes) == head1);
  
  for(e = EdgetreeMinimum(nwp->outedges, head2);
      (tail2 = nwp->outedges[e].value) != 0; /* loop if */
      e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of head2 */
    {
      edges2[nedge2] = tail2;
      ++nedge2;
    }
  for(e = EdgetreeMinimum(nwp->inedges, head2);
      (tail2 = nwp->inedges[e].value) != 0; /* loop if */
      e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of head2 */
    {
      edges2[nedge2] = tail2;
      ++nedge2;
    }
  
  ntoggles = 0;
  for(k=0; k < nedge1; k++){
    if (head1 > edges1[k])
      {
	MHp->togglehead[ntoggles] = edges1[k];
	MHp->toggletail[ntoggles] = head1;
      }
    if (head1 < edges1[k]){
      MHp->togglehead[ntoggles] = head1;
      MHp->toggletail[ntoggles] = edges1[k];
    }
    if(head1 != edges1[k]) ntoggles++;
  }
  
  for(k=0; k < nedge2; k++){
    if (head1 > edges2[k])
      {
	MHp->togglehead[ntoggles] = edges2[k];
	MHp->toggletail[ntoggles] = head1;
      }
    if (head1 < edges2[k]){
      MHp->togglehead[ntoggles] = head1;
      MHp->toggletail[ntoggles] = edges2[k];
    }
    if(head1 != edges2[k]) ntoggles++;
  }
  
  for(k=0; k < nedge2; k++){
    if (head2 > edges2[k])
      {
	MHp->togglehead[ntoggles] = edges2[k];
	MHp->toggletail[ntoggles] = head2;
      }
    if (head2 < edges2[k]){
      MHp->togglehead[ntoggles] = head2;
      MHp->toggletail[ntoggles] = edges2[k];
    }
    if(head2 != edges2[k]) ntoggles++;
  }
  
  for(k=0; k < nedge1; k++){
    if (head2 > edges1[k])
      {
	MHp->togglehead[ntoggles] = edges1[k];
	MHp->toggletail[ntoggles] = head2;
      }
    if (head2 < edges1[k]){
      MHp->togglehead[ntoggles] = head2;
      MHp->toggletail[ntoggles] = edges1[k];
    }
    if(head2 != edges1[k]) ntoggles++;
  }
  free(edges1);
  free(edges2);
}

/*********************
 void MH_ConstantEdgesToggles
*********************/
void MH_ConstantEdgesToggles (MHproposal *MHp, DegreeBound *bd, Network *nwp)  {  
  int noutedge=0, ninedge=0, k, fvalid=0;
  int k0, k1;
  Vertex e, alter, head, tail, tail1;
  MHp->ratio=1.0;
  
  while(noutedge+ninedge==0){
    /* select a node at random */
    head = 1 + unif_rand() * nwp->nnodes;
    ninedge  = nwp->indegree[head];
    noutedge = nwp->outdegree[head];
  }
  
  k0 = (int)(unif_rand() * (noutedge+ninedge)); 
  if (k0 < noutedge){
    k=0;
    for(e = EdgetreeMinimum(nwp->outedges, head);
	((tail = nwp->outedges[e].value) != 0 && k<k0);
	e = EdgetreeSuccessor(nwp->outedges, e)){++k;}
  }else{
    k=0;
    for(e = EdgetreeMinimum(nwp->inedges, head);
	((tail = nwp->inedges[e].value) != 0 && k<(k0-noutedge));
	e = EdgetreeSuccessor(nwp->inedges, e)){++k;}
  }
  
  if ( (!nwp->directed_flag && head > tail) ||
       (nwp->directed_flag && k0 >= noutedge) )
    {
      MHp->togglehead[0] = tail;
      MHp->toggletail[0] = head;
    }else{
      MHp->togglehead[0] = head;
      MHp->toggletail[0] = tail;
    }
  
  k1=0;
  fvalid=0;
  while(fvalid==0 && k1 < 100){
    while((alter = 1 + unif_rand() * nwp->nnodes) == head);
    fvalid=1;
    if(alter == tail){fvalid=0;}
    if (k0 < noutedge || !(nwp->directed_flag)){
      for(e = EdgetreeMinimum(nwp->outedges, head);
	  (fvalid==1 && ((tail1 = nwp->outedges[e].value) != 0));
	  e = EdgetreeSuccessor(nwp->outedges, e)){
	if(alter==tail1){fvalid=0;}}
    }
    if (k0 >= noutedge || !(nwp->directed_flag)){
      for(e = EdgetreeMinimum(nwp->inedges, head);
	  (fvalid==1 && ((tail1 = nwp->inedges[e].value) != 0));
	  e = EdgetreeSuccessor(nwp->inedges, e)){
	if(alter==tail1){fvalid=0;}}
    }
    k1++;
  }
  if (k1 == 100){
    MHp->togglehead[0] = MHp->toggletail[0] = 0;
    MHp->togglehead[1] = MHp->toggletail[1] = 0;
  }
  
  if ( (!nwp->directed_flag && alter > head) ||
       (nwp->directed_flag && k0 < noutedge) )
    {
      MHp->togglehead[1] = head;
      MHp->toggletail[1] = alter;
    }else{
      MHp->togglehead[1] = alter;
      MHp->toggletail[1] = head;
    }
  
  if (!fvalid){
    MHp->togglehead[0] = MHp->toggletail[0] = 0;
    MHp->togglehead[1] = MHp->toggletail[1] = 0;
  }else{  
  }
}

/*********************
 void MH_CondDegSwitchToggles
*********************/
void MH_CondDegSwitchToggles (MHproposal *MHp, DegreeBound *bd, Network *nwp)  {  
  int noutedge, ninedge, i;
  int k, k0, toomany;
  Vertex e, head, tail;
  MHp->ratio=1.0;
  
  /* select a node at random */
  for (i = 0; i < 2; i++){
    toomany=0;
    noutedge=0;
    ninedge=0;
    while(noutedge==0 && ninedge==0 && toomany < 100){
      head = 1 + unif_rand() * nwp->nnodes;
      ninedge=0;
      noutedge=0;
      while(noutedge+ninedge==0){
	/* select a node at random */
	head = 1 + unif_rand() * nwp->nnodes;
	ninedge = nwp->indegree[head];
	noutedge = nwp->outdegree[head];
      }
      ++toomany;
    }
    
    if (toomany == 100){
      MHp->togglehead[0] = MHp->toggletail[0] = 0;
      MHp->togglehead[1] = MHp->toggletail[1] = 0;
    }
    
    k0 = (int)(unif_rand() * (noutedge+ninedge)); 
    if (k0 < noutedge){
      k=0;
      for(e = EdgetreeMinimum(nwp->outedges, head);
	  ((tail = nwp->outedges[e].value) != 0 && k<k0);
	  e = EdgetreeSuccessor(nwp->outedges, e)){++k;}
    }else{
      k=0;
      for(e = EdgetreeMinimum(nwp->inedges, head);
	  ((tail = nwp->inedges[e].value) != 0 && k<(k0-noutedge));
	  e = EdgetreeSuccessor(nwp->inedges, e)){++k;}
    }
    if ( (!nwp->directed_flag && head > tail) ||
	 (nwp->directed_flag && k0 >= noutedge) )
      {
	MHp->togglehead[i] = tail;
	MHp->toggletail[i] = head;
      }else{
	MHp->togglehead[i] = head;
	MHp->toggletail[i] = tail;
      }
  }
  
  if (EdgetreeSearch( MHp->togglehead[0],MHp->toggletail[1], nwp->outedges) ||
      EdgetreeSearch( MHp->togglehead[1],MHp->toggletail[0], nwp->outedges) ){
    MHp->togglehead[0] = MHp->toggletail[0] = 0;
    MHp->togglehead[1] = MHp->toggletail[1] = 0;
  }
  
  if ( (!nwp->directed_flag && MHp->togglehead[0] > MHp->toggletail[1]) )
    {
      MHp->togglehead[2] = MHp->toggletail[1];
      MHp->toggletail[2] = MHp->togglehead[0];
    }else{
      MHp->togglehead[2] = MHp->togglehead[0];
      MHp->toggletail[2] = MHp->toggletail[1];
    }
  
  if ( (!nwp->directed_flag && MHp->togglehead[1] > MHp->toggletail[0]) )
    {
      MHp->togglehead[3] = MHp->toggletail[0];
      MHp->toggletail[3] = MHp->togglehead[1];
    }else{
      MHp->togglehead[3] = MHp->togglehead[1];
      MHp->toggletail[3] = MHp->toggletail[0];
    }
}


