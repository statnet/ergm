/*  File src/MHproposals.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2019 Statnet Commons
 */
#include "MHproposals.h"
#include "ergm_edgelist.h"
#include "ergm_changestat.h"
#include "ergm_rlebdm.h"
#include "ergm_MHstorage.h"

/*********************
 void MH_randomtoggle

 Default MH algorithm
*********************/
MH_P_FN(MH_randomtoggle){  

  /* *** don't forget tail-> head now */

  if(MHp->ntoggles == 0) { /* Initialize randomtoggle */
    MHp->ntoggles=1;
    return;
  }
  
  BD_LOOP({
      GetRandDyad(Mtail, Mhead, nwp);
    });
}

/********************
   void MH_TNT
   Tie/no tie:  Gives at least 50% chance of
   proposing a toggle of an existing edge, as opposed
   to simple random toggles that rarely do so in sparse 
   networks
***********************/
MH_P_FN(MH_TNT)
{
  /* *** don't forget tail-> head now */
  
  Edge nedges=EDGECOUNT(nwp);
  static double comp=0.5;
  static double odds;
  static Dyad ndyads;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    odds = comp/(1.0-comp);
    ndyads = DYADCOUNT(nwp);
    return;
  }

  double logratio=0;
  BD_LOOP({
      if (unif_rand() < comp && nedges > 0) { /* Select a tie at random */
	GetRandEdge(Mtail, Mhead, nwp);
	/* Thanks to Robert Goudie for pointing out an error in the previous 
	   version of this sampler when proposing to go from nedges==0 to nedges==1 
	   or vice versa.  Note that this happens extremely rarely unless the 
	   network is small or the parameter values lead to extremely sparse 
	   networks.  */
	logratio = log((nedges==1 ? 1.0/(comp*ndyads + (1.0-comp)) :
			      nedges / (odds*ndyads + nedges)));
      }else{ /* Select a dyad at random */
	GetRandDyad(Mtail, Mhead, nwp);
	
	if(IS_OUTEDGE(Mtail[0],Mhead[0])!=0){
	  logratio = log((nedges==1 ? 1.0/(comp*ndyads + (1.0-comp)) :
				nedges / (odds*ndyads + nedges)));
	}else{
	  logratio = log((nedges==0 ? comp*ndyads + (1.0-comp) :
				1.0 + (odds*ndyads)/(nedges + 1)));
	}
      }
    });
  MHp->logratio += logratio;
}

/********************
    MH_StratTNT
********************/

MH_I_FN(Mi_StratTNT) {
  // process the inputs and initialize all the edgelists in storage; set MHp->ntoggles to 1
  MHp->ntoggles = 1;
  
  int nmixtypes = MHp->inputs[0];
  
  double *tailtypes = MHp->inputs + 1;
  double *headtypes = tailtypes + nmixtypes;
  
  int nattrcodes = MHp->inputs[1 + 3*nmixtypes];
  
  double *vattr = MHp->inputs + 1 + 3*nmixtypes + 1 + nattrcodes + N_NODES;
  
  // in storage, we need:
  // (1,2) integer count of how many edges each mixing type has, and what the maximum allowed in current storage is
  // (3,4) tail and head vertex lists for edges of each mixing type
  // (5,6) when adding/removing an edge in the U function below, it will save time to know its mixing type and 
  //       (if applicable) index within the edgelist so the last two entries will be an int * mixingtype and an int * edgeindex
  ALLOC_STORAGE(6, void *, sto);
  
  Vertex **tailsvec = (Vertex **)Calloc(nmixtypes, Vertex *);
  Vertex **headsvec = (Vertex **)Calloc(nmixtypes, Vertex *);
  int *mixcounts = (int *)Calloc(nmixtypes, int);
  int *maxcounts = (int *)Calloc(nmixtypes, int);
    
  for(int i = 0; i < nmixtypes; i++) {
    // start with maxcounts[i] = nwp->maxedges, which is assumed to be larger than the number 
    // of edges in the network at the time this initialization function is called
    tailsvec[i] = (Vertex *)Calloc(nwp->maxedges, Vertex);
    headsvec[i] = (Vertex *)Calloc(nwp->maxedges, Vertex);
    maxcounts[i] = nwp->maxedges;
  }
  
  Vertex head;
  Edge e;
  for(Vertex tail = 1; tail <= N_NODES; tail++) {
    STEP_THROUGH_OUTEDGES(tail, e, head) {
      // i indexes possible mixing types
      for(int i = 0; i < nmixtypes; i++) {
        // in undirected unipartite case we need to allow the combination of attr values to occur in either order; in other cases (dir. or bip.), mixing is directed
        if((vattr[tail - 1] == tailtypes[i] && vattr[head - 1] == headtypes[i]) || (!DIRECTED && !BIPARTITE && vattr[head - 1] == tailtypes[i] && vattr[tail - 1] == headtypes[i])) {
          // this edge is of mixing type i and should be added to the ith edgelist, updating the count accordingly
          // note that tail < head is guaranteed in the undirected case since e is an outedge from tail to head
          tailsvec[i][mixcounts[i]] = tail;
          headsvec[i][mixcounts[i]] = head;
          mixcounts[i]++;
          break; // so we don't bother checking remaining mixing types for this particular edge
        }
      }
    }
  }
  
  // assign counts and edgelists to storage
  sto[0] = (void *)mixcounts;
  sto[1] = (void *)maxcounts;
  
  sto[2] = (void *)tailsvec;
  sto[3] = (void *)headsvec;

  // make room for mixingtype and edgeindex
  sto[4] = Calloc(1, int);
  sto[5] = Calloc(1, int);  
}

MH_P_FN(MH_StratTNT) {
  int nmixtypes = MHp->inputs[0];

  double *pmat = MHp->inputs + 1 + 2*nmixtypes;
  
  double ur = unif_rand();
  
  // find the first mixing type i with (cumulative) probability larger than ur
  int i = 0;
  while(ur > pmat[i]) {
    i++;
  }
  
  GET_STORAGE(void *, sto);

  // record the mixing type of the toggle, in case it's needed in the U function later
  ((int *)sto[4])[0] = i;    
  
  int *mixcounts = (int *)sto[0];
  
  Vertex **tailsvec = (Vertex **)sto[2];
  Vertex **headsvec = (Vertex **)sto[3];

  // for the logratio calculation, we will need some statistics
  double *tailtypes = MHp->inputs + 1;
  double *headtypes = tailtypes + nmixtypes;
  
  int tailtype = tailtypes[i];
  int headtype = headtypes[i];
  
  double *nodecountsbycode = MHp->inputs + 1 + 3*nmixtypes + 1;  
  
  // number of edges of this mixing type
  int nedgestype = mixcounts[i];

  // number of dyads of this mixing type
  // note tailtype == headtype will never hold
  // for bipartite graphs, given the way
  // we recoded things in R
  int ndyadstype;
  if(tailtype == headtype) {
    if(DIRECTED)
      ndyadstype = nodecountsbycode[tailtype - 1]*(nodecountsbycode[headtype - 1] - 1);
    else
      ndyadstype = nodecountsbycode[tailtype - 1]*(nodecountsbycode[headtype - 1] - 1)/2;
  } else {
    ndyadstype = nodecountsbycode[tailtype - 1]*nodecountsbycode[headtype - 1];
  }
  
  double logratio = 0;

  BD_LOOP({
    if(unif_rand() < 0.5 && mixcounts[i] > 0) {
      // select an existing edge of type i at random, and propose toggling it off
      int edgeindex = mixcounts[i]*unif_rand();
      
      // in the undirected case, we always have tailsvec[i][edgeindex] < headsvec[i][edgeindex],
      // so we don't need to check that here
      Mtail[0] = tailsvec[i][edgeindex];
      Mhead[0] = headsvec[i][edgeindex];
      
      // set edgeindex in storage in case we need to update
      ((int *)sto[5])[0] = edgeindex;
      
      // logratio is essentially copied from TNT, because the probability of 
      // choosing this particular mixing type cancels upon taking the ratio;
      // still need to count only edges and dyads of the appropriate mixing type, though
  	  logratio = log((nedgestype == 1 ? 1.0/(0.5*ndyadstype + 0.5) :
                      nedgestype / ((double) ndyadstype + nedgestype)));
    } else {
      // select a dyad of type i and propose toggling it
      int nattrcodes = MHp->inputs[1 + 3*nmixtypes];
      
      double *nodeindicesbycode = MHp->inputs + 1 + 3*nmixtypes + 1 + nattrcodes;
  
      int tailindex = nodecountsbycode[tailtype - 1]*unif_rand();
      int headindex;
      if(tailtype == headtype) {
        // need to avoid sampling a loop
        headindex = (nodecountsbycode[headtype - 1] - 1)*unif_rand();
        if(headindex == tailindex) {
          headindex = nodecountsbycode[headtype - 1] - 1;
        }
      } else {
        // any old head will do
        headindex = nodecountsbycode[headtype - 1]*unif_rand();
      }
            
      for(int j = 0; j < tailtype - 1; j++) {
        tailindex += nodecountsbycode[j];
      }
        
      for(int j = 0; j < headtype - 1; j++) {
        headindex += nodecountsbycode[j];
      }
      
      Vertex tail = nodeindicesbycode[tailindex];
      Vertex head = nodeindicesbycode[headindex];
      
      if(tail > head && !DIRECTED) {
        Vertex tmp = tail;
        tail = head;
        head = tmp;
      }
      
      if(IS_OUTEDGE(tail,head)) {
        // pick a new edge from the edgelist uniformly at random so we know its index
        // and hence don't have to look up the index of the edge tail -> head; this gives
        // the same probability of picking each existing edge as if we used the tail -> head
        // edge, but also allows us to keep the edgelists unsorted (at the cost of generating
        // an extra random index in this case)
        int edgeindex = unif_rand()*mixcounts[i];
                
        Mtail[0] = tailsvec[i][edgeindex];
        Mhead[0] = headsvec[i][edgeindex];

        ((int *)sto[5])[0] = edgeindex;

        logratio = log((nedgestype == 1 ? 1.0/(0.5*ndyadstype + 0.5) :
                        nedgestype / ((double) ndyadstype + nedgestype)));
      }else{
        Mtail[0] = tail;
        Mhead[0] = head;
        
        ((int *)sto[5])[0] = -1; // meaning not currently in edgelist
                
        logratio = log((nedgestype == 0 ? 0.5*ndyadstype + 0.5 :
                        1.0 + (ndyadstype)/((double) nedgestype + 1)));
      }
    }
  });
  
  MHp->logratio += logratio;
}

MH_U_FN(Mu_StratTNT) {
  // add or remove edge from appropriate edgelist; update edgelist length
  GET_STORAGE(void *, sto);

  int *mixcounts = (int *)sto[0];
  int *maxcounts = (int *)sto[1];
  
  Vertex **tailsvec = (Vertex **)sto[2];
  Vertex **headsvec = (Vertex **)sto[3];

  int mixtype = ((int *)sto[4])[0];
  int edgeindex = ((int *)sto[5])[0];
  
  if(edgeindex >= 0) {
    // we are removing an existing edge; move the last edge into the place of the edge we're removing;
    // this movement is unnecessary if we are removing the last edge, but oh well (it doesn't hurt anything)
    tailsvec[mixtype][edgeindex] = tailsvec[mixtype][mixcounts[mixtype] - 1];
    headsvec[mixtype][edgeindex] = headsvec[mixtype][mixcounts[mixtype] - 1];

    // decrement the count
    mixcounts[mixtype]--;
  } else {
    // we are adding a new edge; put it on the end of the edgelist
    
    // check if mixcounts[mixtype] is maxcounts[mixtype];
    // if so, copy everything to a new location twice the size, and *then* add the new edge
    // being sure to update storage pters appropriately    
    if(mixcounts[mixtype] == maxcounts[mixtype]) {
      Vertex *newtailsvec = (Vertex *)Calloc(2*maxcounts[mixtype], Vertex);
      Vertex *newheadsvec = (Vertex *)Calloc(2*maxcounts[mixtype], Vertex);
      memcpy(newtailsvec, tailsvec[mixtype], sizeof(Vertex)*maxcounts[mixtype]);
      memcpy(newheadsvec, headsvec[mixtype], sizeof(Vertex)*maxcounts[mixtype]);
      maxcounts[mixtype] = 2*maxcounts[mixtype];
      Free(tailsvec[mixtype]);
      Free(headsvec[mixtype]);
      tailsvec[mixtype] = newtailsvec;
      headsvec[mixtype] = newheadsvec;
    }

    tailsvec[mixtype][mixcounts[mixtype]] = tail;
    headsvec[mixtype][mixcounts[mixtype]] = head;
    
    // increment the count
    mixcounts[mixtype]++;
  }
}

MH_F_FN(Mf_StratTNT) {
  // Free all the things
  GET_STORAGE(void *, sto);
  
  int nmixtypes = MHp->inputs[0];
    
  Free(sto[5]);
  Free(sto[4]);

  Vertex **headsvec = (Vertex **)sto[3];  
  Vertex **tailsvec = (Vertex **)sto[2];

  for(int i = 0; i < nmixtypes; i++) {
    Free(headsvec[i]);
    Free(tailsvec[i]);
  }

  Free(headsvec);  
  Free(tailsvec);

  Free(sto[1]);  
  Free(sto[0]);
  
  // MHp->storage itself should be Freed by MHProposalDestroy
}

/********************
   void MH_TNT10
   Attempts to do 10 TNT steps at once, but this seems flawed currently
   because it does not correctly update network quantities like nedges
   after each of the 10 proposed toggles.
***********************/
MH_P_FN(MH_TNT10)
{
  /* *** don't forget tail-> head now */
  
  Edge nedges=EDGECOUNT(nwp);
  static double comp=0.5;
  static double odds;
  static Dyad ndyads;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=10;
    odds = comp/(1.0-comp);
    ndyads = DYADCOUNT(nwp);
    return;
  }
  
  double logratio = 0;
  BD_LOOP({
      logratio = 0;
      for(unsigned int n = 0; n < 10; n++){
	if (unif_rand() < comp && nedges > 0) { /* Select a tie at random */
	  GetRandEdge(Mtail, Mhead, nwp);
	  logratio += log(nedges  / (odds*ndyads + nedges));
	}else{ /* Select a dyad at random */
	  GetRandDyad(Mtail+n, Mhead+n, nwp);
	  if(IS_OUTEDGE(Mtail[n],Mhead[n])!=0){
	    logratio += log((nedges==1 ? 1.0/(comp*ndyads + (1.0-comp)) :
				  nedges / (odds*ndyads + nedges)));
	  }else{
	    logratio += log((nedges==0 ? comp*ndyads + (1.0-comp) :
				  1.0 + (odds*ndyads)/(nedges + 1)));
	  }
	} 
      }
    });
  MHp->logratio += logratio;
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
MH_P_FN(MH_ConstantEdges){  
  /* *** don't forget tail-> head now */
  
  if(MHp->ntoggles == 0) { /* Initialize */
    if(nwp->nedges==0 || nwp->nedges==DYADCOUNT(nwp)) MHp->ntoggles=MH_FAILED; /* Empty or full network. */
    else MHp->ntoggles=2;
    return;
  }
  /* Note:  This proposal cannot be used for full or empty observed graphs.
     If desired, we could check for this at initialization phase. 
     (For now, however, no way to easily return an error message and stop.)*/
  BD_LOOP({
      /* First, select edge at random */
      GetRandEdge(Mtail, Mhead, nwp);
      /* Second, select non-edge at random */
      GetRandNonedge(Mtail+1, Mhead+1, nwp);
    });
}

/*********************
 void MH_CondDegreeDist
 It used to be called  MH_CondDegDistSwapToggles
*********************/
MH_P_FN(MH_CondDegreeDist){  
  int noutedge=0, ninedge=0, k, fvalid;
  int k0, j0, j1, k1;
  int j0h, j1h;
  int trynode;
  Vertex e, alter, tail=0, head, head1;
  
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
    tail = 1 + unif_rand() * N_NODES;
    ninedge  = nwp->indegree[tail];
    noutedge = nwp->outdegree[tail];
  }

  /* choose a edge of the node at random */
    /* *** don't forget tail-> head now */

  k0 = (int)(unif_rand() * (noutedge+ninedge)); 
  if (k0 < noutedge){
    k=0;
    for(e = EdgetreeMinimum(nwp->outedges, tail);
    ((head = nwp->outedges[e].value) != 0 && k<k0);
    e = EdgetreeSuccessor(nwp->outedges, e)){++k;}
  }else{
    k=0;
    for(e = EdgetreeMinimum(nwp->inedges, tail);
    ((head = nwp->inedges[e].value) != 0 && k<(k0-noutedge));
    e = EdgetreeSuccessor(nwp->inedges, e)){++k;}
  }

  if ( (!DIRECTED && tail > head) ||
  (DIRECTED && k0 >= noutedge) ) {
    Mtail[0] = head;
    Mhead[0] = tail;
  }else{
    Mtail[0] = tail;
    Mhead[0] = head;
  }
  
  k1=0;
  fvalid=0;
  while(fvalid==0 && k1 < 100){
    while((alter = 1 + unif_rand() * N_NODES) == tail);
    fvalid=1;
    if(alter == head){fvalid=0;}
    if (k0 < noutedge || !DIRECTED){
      for(e = EdgetreeMinimum(nwp->outedges, tail);
      (fvalid==1 && ((head1 = nwp->outedges[e].value) != 0));
      e = EdgetreeSuccessor(nwp->outedges, e)){
        if(alter==head1){fvalid=0;}}
    }
    if (k0 >= noutedge || !DIRECTED){
      for(e = EdgetreeMinimum(nwp->inedges, tail);
      (fvalid==1 && ((head1 = nwp->inedges[e].value) != 0));
      e = EdgetreeSuccessor(nwp->inedges, e)){
        if(alter==head1){fvalid=0;}}
    }
    k1++;
  }

  if (k1 == 100){
    fvalid=0;
    continue;
  }
  
  if ( (!DIRECTED && alter > tail) ||
       (DIRECTED && k0 < noutedge) )
    {
      Mtail[1] = tail;
      Mhead[1] = alter;
    }else{
      Mtail[1] = alter;
      Mhead[1] = tail;
    }
  
  if (!DIRECTED){
    /* Check undirected degrees */
    k0 =nwp->outdegree[tail]  + nwp->indegree[tail];
    j0h=nwp->outdegree[head]  + nwp->indegree[head];
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
     j0h=nwp->indegree[head];
     j1h=nwp->indegree[alter];
   }else{
     /* Check outdegrees */
     j0h=nwp->outdegree[head];
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
    Mtail[1] = Mtail[0];
    Mhead[1] = Mhead[0];
  }
}

/*********************
 void MH_CondOutDegreeDist
*********************/
MH_P_FN(MH_CondOutDegreeDist){  
  int noutedge=0, k, fvalid=0;
  int k0, k1;
  int trynode;
  Vertex e, alter, tail=0, head, head1;
  
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
    tail = 1 + unif_rand() * N_NODES;
    noutedge = nwp->outdegree[tail];
  }
  
  k0 = (int)(unif_rand() * noutedge); 
  k=0;
  for(e = EdgetreeMinimum(nwp->outedges, tail);
      ((head = nwp->outedges[e].value) != 0 && k<k0);
      e = EdgetreeSuccessor(nwp->outedges, e)){++k;}
  Mtail[0] = tail;
  Mhead[0] = head;
  
  k1=0;
  fvalid=0;
  while(fvalid==0 && k1 < 100){
    while((alter = 1 + unif_rand() * N_NODES) == tail);
    fvalid=1;
    if(alter == head){fvalid=0;}
    for(e = EdgetreeMinimum(nwp->outedges, tail);
	(fvalid==1 && ((head1 = nwp->outedges[e].value) != 0));
	e = EdgetreeSuccessor(nwp->outedges, e)){
      if(alter==head1){fvalid=0;}}
    k1++;
  }
  if (k1 == 100){
    fvalid=0;
    continue;
  }
  
  Mtail[1] = tail;
  Mhead[1] = alter;
  }
  
  if(trynode==1500 || !CheckTogglesValid(MHp, nwp)){
      Mtail[0] = 1;
      Mhead[0] = 2;
      Mtail[1] = 1;
      Mhead[1] = 2;
  }
  

}

/*********************
 void MH_CondInDegreeDist
*********************/
MH_P_FN(MH_CondInDegreeDist){  
  int ninedge=0, k, fvalid=0;
  int k0, k1;
  int trynode;
  Vertex e, alter, tail=0, head, head1;

  /* *** don't forget tail-> head now */

  
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
    tail = 1 + unif_rand() * N_NODES;
    ninedge = nwp->indegree[tail];
  }
  
  k0 = (int)(unif_rand() * ninedge); 
  k=0;
  for(e = EdgetreeMinimum(nwp->inedges, tail);
      ((head = nwp->inedges[e].value) != 0 && k<k0);
      e = EdgetreeSuccessor(nwp->inedges, e)){++k;}
  Mtail[0] = head;
  Mhead[0] = tail;
  
  k1=0;
  fvalid=0;
  while(fvalid==0 && k1 < 100){
    while((alter = 1 + unif_rand() * N_NODES) == tail);
    fvalid=1;
    if(alter == head){fvalid=0;}
    for(e = EdgetreeMinimum(nwp->inedges, tail);
	(fvalid==1 && ((head1 = nwp->inedges[e].value) != 0));
	e = EdgetreeSuccessor(nwp->inedges, e)){
      if(alter==head1){fvalid=0;}}
    k1++;
  }
  if (k1 == 100){
    fvalid=0;
    continue;
  }
  
  Mtail[1] = alter;
  Mhead[1] = tail;
  
  }
  
  if(trynode==1500){
      Mtail[0] = 1;
      Mhead[0] = 2;
      Mtail[1] = 1;
      Mhead[1] = 2;
  }
}

/*********************
 void MH_TwoRandomToggles
*********************/
MH_P_FN(MH_TwoRandomToggles){  
  Vertex tail, head;
  int i;

  /* *** don't forget tail-> head now */
  
  if(MHp->ntoggles == 0) { /* Initialize OneRandomToggle */
    MHp->ntoggles=2;
    return;
  }

  for (i = 0; i < 2; i++){
   tail = 1 + unif_rand() * N_NODES;
   while ((head = 1 + unif_rand() * N_NODES) == tail);
   if (!DIRECTED && tail > head) {
     Mtail[i] = head;
     Mhead[i] = tail;
   }else{
     Mtail[i] = tail;
     Mhead[i] = head;
   }
  }
}

/*********************
 void MH_RandomNode
*********************/
MH_P_FN(MH_randomnode){
  
  Vertex root, alter;
  int j;
  
  if(MHp->ntoggles == 0) { /* Initialize OneRandomToggle */
    MHp->ntoggles= N_NODES - 1;
    return;
  }

  root = 1 + unif_rand() * N_NODES;
  
  j = 0;
  for (alter = 1; alter <= N_NODES; alter++)
    {
      /* there is never an edge (root, root) */
      if (alter != root) {
       if (!DIRECTED && root > alter) {
        Mtail[j] = alter;
        Mhead[j] = root;
       }else{
        Mtail[j] = root;
        Mhead[j] = alter;
       }
       j++;
      }
    }
}

/********************
   void MH_randomtoggleList
   Propose ONLY edges on a static list
***********************/
MH_P_FN(MH_randomtoggleList)
{  
  static Edge nedges0;

  if(MHp->ntoggles == 0) { /* Initialize */
    nedges0 = MH_INPUTS[0];
    if(nedges0==0) MHp->ntoggles=MH_FAILED; /* Dyad list has no elements. */
    else MHp->ntoggles=1;
    return;
  }
  
  BD_LOOP({
      /* Select a dyad at random that is in the reference graph. (We
	 have a convenient sampling frame.) */
      /* Generate. */
      Edge rane = 1 + unif_rand() * nedges0;
      Mtail[0]=MH_INPUTS[rane];
      Mhead[0]=MH_INPUTS[nedges0+rane];
    });
}

/********************
   void MH_randomtoggleRLE
   Propose ONLY edges on an RLE-compressed list
***********************/
MH_I_FN(Mi_RLE){
  ALLOC_STORAGE(1, RLEBDM1D, r);
  double *inputs = MHp->inputs;
  *r = unpack_RLEBDM1D(&inputs, nwp->nnodes);
  if(r->ndyads==0) MHp->ntoggles=MH_FAILED; /* Dyad list has no elements. */
  else MHp->ntoggles=1;
}

MH_P_FN(Mp_RLE){
  GET_STORAGE(RLEBDM1D, r);

  BD_LOOP({
      /* Select a dyad at random that is in the reference graph. (We
	 have a convenient sampling frame.) */
      /* Generate. */
      GetRandRLEBDM1D_RS(Mtail, Mhead, r);
    });
}

/********************
   void MH_listTNT
   Propose ONLY edges on a static list
   Use TNT weights.
   This is a fusion of MH_DissolutionMLETNT and MH_TNT:

   The "intersect" network is requested that is the intersection of
   dyads on the static list and the edges present in nwp. Then,
   standard TNT procedure is followed, but the dyad space (and the
   number of dyads) is the number of dyads in the static list and the
   network for the ties is the ties in the intersect network.
***********************/
MH_I_FN(Mi_listTNT){
  Dyad ndyads = MH_INPUTS[0]; // Note that ndyads here is the number of dyads in the list.
  if(ndyads==0){
    MHp->ntoggles=MH_FAILED; /* Dyad list has no elements. */
    return;
  }else MHp->ntoggles=1;
  double *list = MH_INPUTS+1;
  Network *intersect = STORAGE = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, 0, 0, NULL);
  for(Edge i=0; i<ndyads; i++){
    Vertex tail=list[i], head=list[ndyads+i];
    if(IS_OUTEDGE(tail, head)!=0)
      ToggleEdge(tail, head, intersect);
  }
}

MH_U_FN(Mu_listTNT){
  Network *intersect = STORAGE;
  ToggleEdge(tail, head, intersect);
}

MH_P_FN(Mp_listTNT){
  const double comp=0.5, odds=comp/(1.0-comp);
  Dyad ndyads = MH_INPUTS[0]; // Note that ndyads here is the number of dyads in the list.
  double *list = MH_INPUTS+1;

  Network *intersect = STORAGE;

  Edge nedges=EDGECOUNT(intersect);
  
  double logratio=0;
  BD_LOOP({
      if (unif_rand() < comp && nedges > 0) { /* Select a tie at random from the network of eligibles */
	GetRandEdge(Mtail, Mhead, intersect);
	/* Thanks to Robert Goudie for pointing out an error in the previous 
	   version of this sampler when proposing to go from nedges==0 to nedges==1 
	   or vice versa.  Note that this happens extremely rarely unless the 
	   network is small or the parameter values lead to extremely sparse 
	   networks.  */
	logratio = log((nedges==1 ? 1.0/(comp*ndyads + (1.0-comp)) :
			      nedges / (odds*ndyads + nedges)));
      }else{ /* Select a dyad at random from the list */
	Edge rane = unif_rand() * ndyads;
	Mtail[0]=list[rane];
	Mhead[0]=list[ndyads+rane];
	
	if(IS_OUTEDGE(Mtail[0],Mhead[0],intersect)){
	  logratio = log((nedges==1 ? 1.0/(comp*ndyads + (1.0-comp)) :
				nedges / (odds*ndyads + nedges)));
	}else{
	  logratio = log((nedges==0 ? comp*ndyads + (1.0-comp) :
				1.0 + (odds*ndyads)/(nedges + 1)));
	}
      }
    });
  MHp->logratio += logratio;
}

MH_F_FN(Mf_listTNT){
  Network *intersect = STORAGE;
  NetworkDestroy(intersect);
  STORAGE = NULL;
}

/********************
   void MH_RLETNT
   Propose ONLY edges on a static list
   Use TNT weights.
   This is a fusion of MH_DissolutionMLETNT and MH_TNT:

   A "intersect" network is constructed that is the intersection of
   dyads on the static list and the edges present in nwp. Then,
   standard TNT procedure is followed, but the dyad space (and the
   number of dyads) is the number of dyads in the static list and the
   network for the ties is the ties in the discord network.
***********************/
typedef struct {
  RLEBDM1D r;
  Network *intersect;
} StoreRLEBDM1DAndNet;

MH_I_FN(Mi_RLETNT){
  ALLOC_STORAGE(1, StoreRLEBDM1DAndNet, storage);
  double *inputs = MHp->inputs;
  storage->r = unpack_RLEBDM1D(&inputs, nwp->nnodes);
  if(storage->r.ndyads==0){
    MHp->ntoggles=MH_FAILED; /* Dyad list has no elements. */
    return;
  }else MHp->ntoggles=1;
  storage->intersect = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, 0, 0, NULL);
  EXEC_THROUGH_NET_EDGES_PRE(t, h, e, {
      if(GetRLEBDM1D(t, h, &storage->r)){
	ToggleEdge(t, h, storage->intersect);
      }
    });
  
  if(EDGECOUNT(storage->intersect)==EDGECOUNT(nwp)){ // There are no ties in the initial network that are fixed.
    NetworkDestroy(storage->intersect);
    storage->intersect = NULL; // "Signal" that there is no discordance network.
  }
}

MH_P_FN(Mp_RLETNT){
  GET_STORAGE(StoreRLEBDM1DAndNet, storage);

  const double comp=0.5, odds=comp/(1.0-comp);

  Network *nwp1 = storage->intersect ? storage->intersect : nwp;
  Edge nedges= EDGECOUNT(nwp1);
  double logratio=0;
  BD_LOOP({
      if (unif_rand() < comp && nedges > 0) { /* Select a tie at random from the network of eligibles */
	GetRandEdge(Mtail, Mhead, nwp1);
	/* Thanks to Robert Goudie for pointing out an error in the previous 
	   version of this sampler when proposing to go from nedges==0 to nedges==1 
	   or vice versa.  Note that this happens extremely rarely unless the 
	   network is small or the parameter values lead to extremely sparse 
	   networks.  */
	logratio = log((nedges==1 ? 1.0/(comp*storage->r.ndyads + (1.0-comp)) :
			      nedges / (odds*storage->r.ndyads + nedges)));
      }else{ /* Select a dyad at random from the list */
	GetRandRLEBDM1D_RS(Mtail, Mhead, &storage->r);
	
	if(IS_OUTEDGE(Mtail[0],Mhead[0],nwp1)){
	  logratio = log((nedges==1 ? 1.0/(comp*storage->r.ndyads + (1.0-comp)) :
				nedges / (odds*storage->r.ndyads + nedges)));
	}else{
	  logratio = log((nedges==0 ? comp*storage->r.ndyads + (1.0-comp) :
				1.0 + (odds*storage->r.ndyads)/(nedges + 1)));
	}
      }
    });
  MHp->logratio += logratio;
}

MH_U_FN(Mu_RLETNT){
  GET_STORAGE(StoreRLEBDM1DAndNet, storage);
  if(storage->intersect) ToggleEdge(tail, head, storage->intersect);
}

MH_F_FN(Mf_RLETNT){
  GET_STORAGE(StoreRLEBDM1DAndNet, storage);
  if(storage->intersect) NetworkDestroy(storage->intersect);
}

/* The ones below have not been tested */

/*********************
 void MH_ConstrainedCondOutDegDist
*********************/
MH_P_FN(MH_ConstrainedCondOutDegDist){  
  int noutedge=0, k, fvalid=0;
  int k0, k1;
  Vertex e, alter, tail, head, head1;

  /* *** don't forget tail-> head now */

  while(noutedge==0){
    /* select a node at random */
    tail = 1 + unif_rand() * N_NODES;
    noutedge = nwp->outdegree[tail];
  }
  
  k0 = (int)(unif_rand() * noutedge); 
  k=0;
  for(e = EdgetreeMinimum(nwp->outedges, tail);
      ((head = nwp->outedges[e].value) != 0 && k<k0);
      e = EdgetreeSuccessor(nwp->outedges, e)){++k;}
  Mtail[0] = tail;
  Mhead[0] = head;
  
  k1=0;
  fvalid=0;
  while(fvalid==0 && k1 < 100){
    while((alter = 1 + unif_rand() * N_NODES) == tail);
    fvalid=1;
    if(alter == head){fvalid=0;}
    for(e = EdgetreeMinimum(nwp->outedges, tail);
	(fvalid==1 && ((head1 = nwp->outedges[e].value) != 0));
	e = EdgetreeSuccessor(nwp->outedges, e)){
      if(alter==head1){fvalid=0;}}
    k1++;
  }
  if (k1 == 100){
    Mtail[0] = Mhead[0] = 0;
    Mtail[1] = Mhead[1] = 0;
  }
  
  Mtail[1] = tail;
  Mhead[1] = alter;
  
  if (!fvalid){
    Mtail[0] = Mhead[0] = 0;
    Mtail[1] = Mhead[1] = 0;
  }
  
  for(k=0; k < 2; k++){
    if (dEdgeListSearch(Mtail[k], Mhead[k], MH_INPUTS)==0){
      Mtail[0] = Mhead[0] = 0;
      Mtail[1] = Mhead[1] = 0;
    }
  }
}


MH_P_FN(MH_NodePairedTiesToggles){  
  /* chooses a node and toggles all ties and
	 and toggles an equal number of matching nonties
	 for that node */
  int nedge=0,j,k;
  int fvalid = 1;
  Vertex e, tail, prop;

  /* *** don't forget tail-> head now */
  
  /* double to integer coercion */
  tail = 1 + unif_rand() * N_NODES; 
  
  for(e = EdgetreeMinimum(nwp->outedges, tail);
      (prop = nwp->outedges[e].value) != 0; /* loop if */
      e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
    {
      Mtail[nedge] = tail;
      Mhead[nedge] = prop;
      ++nedge;
    }
  for(e = EdgetreeMinimum(nwp->inedges, tail);
      (prop = nwp->inedges[e].value) != 0; /* loop if */
      e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail */
    {
      Mhead[nedge] = tail;
      Mtail[nedge] = prop;
      ++nedge;
    }
  
  if(nedge > N_NODES-nedge){
    Mtail[0] = Mhead[0] = 0;
    Mtail[1] = Mhead[1] = 0;
  }  
  j = 0;
  while (j <=nedge)
    {
      prop = 1 + unif_rand() * N_NODES; 
      k=0;
      fvalid=1;
      while(fvalid==1 && k<nedge+j){
	if(IS_OUTEDGE( MIN(prop,Mtail[k]),
		       MAX(prop,Mtail[k])) +
	   IS_OUTEDGE( MIN(prop,Mhead[k]),
		       MAX(prop,Mhead[k]))==0
	   ){++k;
	}else{
	  fvalid=0;
	}
      }
      if(prop>tail){
	Mtail[j+nedge] = tail;
	Mhead[j+nedge] = prop;
      }else{
	Mtail[j+nedge] = prop;
	Mhead[j+nedge] = tail;
      }
      ++j;
    }
  
  j = 2*nedge;
  if (!CheckTogglesValid(MHp, nwp))
    {
      *Mtail = *Mhead = 0;
    }
}

/*********************
 void MH_OneRandomTnTNode
*********************/
MH_P_FN(MH_OneRandomTnTNode){  
  Vertex tail=0, head, e, head1;
  int noutedge=0, ninedge=0, k0=0, fvalid=0, k;
  /* int ndyad; */

  /* *** don't forget tail-> head now */
  
  /* if ( DIRECTED )
    {
      ndyad = (N_NODES - 1) * N_NODES;
    }else{
      ndyad = (N_NODES - 1) * N_NODES / 2;
    } */

  double logratio=0;
  fvalid=0;
  while(fvalid==0){
    
    if ( unif_rand() < 0.5 && EDGECOUNT(nwp) > 0) 
      {
	
	/* select a tie */
	ninedge=0;
	noutedge=0;
	while(noutedge+ninedge==0){
	  /* select a node at random */
	  tail = 1 + unif_rand() * N_NODES;
	  ninedge = nwp->indegree[tail];
	  noutedge = nwp->outdegree[tail];
	}
	
	k0 = (int)(unif_rand() * (noutedge+ninedge)); 
	if (k0 < noutedge){
	  k=0;
	  for(e = EdgetreeMinimum(nwp->outedges, tail);
	      ((head = nwp->outedges[e].value) != 0 && k<k0);
	      e = EdgetreeSuccessor(nwp->outedges, e)){++k;}
	}else{
	  k=0;
	  for(e = EdgetreeMinimum(nwp->inedges, tail);
	      ((head = nwp->inedges[e].value) != 0 && k<(k0-noutedge));
	      e = EdgetreeSuccessor(nwp->inedges, e)){++k;}
	}
	if ( (!DIRECTED && tail > head) ||
	     (DIRECTED && k0 >= noutedge) )
	  {
	    Mtail[0] = head;
	    Mhead[0] = tail;
	  }else{
	    Mtail[0] = tail;
	    Mhead[0] = head;
	  }
	
	logratio = log(((noutedge+ninedge)*1.0)/(N_NODES-1-noutedge-ninedge-1));
	fvalid =1;
      }else{
	/* Choose random non-tie */

	/* select a node at random */
	ninedge=N_NODES-1;
	noutedge=0;
	while(noutedge+ninedge>=(N_NODES-1)){
	  ninedge=0;
	  /* select a node at random */
	  tail = 1 + unif_rand() * N_NODES;
	  ninedge = nwp->indegree[tail];
	  noutedge = nwp->outdegree[tail];
	}
	
	fvalid=0;
	while(fvalid==0){
	  while ((head = 1 + unif_rand() * N_NODES) == tail);
	  fvalid=1;
	  for(e = EdgetreeMinimum(nwp->outedges, tail);
	      (fvalid==1 && ((head1 = nwp->outedges[e].value) != 0));
	      e = EdgetreeSuccessor(nwp->outedges, e)){
	    if(head==head1){fvalid=0;}}
	  if (!(DIRECTED)){
	    for(e = EdgetreeMinimum(nwp->inedges, tail);
		(fvalid==1 && ((head1 = nwp->inedges[e].value) != 0));
		e = EdgetreeSuccessor(nwp->inedges, e)){
	      if(head==head1){fvalid=0;}}
	  }
	}
	
	if ( (!DIRECTED && tail > head) ||
	     (DIRECTED && k0 >= noutedge) )
	  {
	    Mtail[0] = head;
	    Mhead[0] = tail;
	  }else{
	    Mtail[0] = tail;
	    Mhead[0] = head;
	  }
	
        if ( DIRECTED )
	  {
	    logratio = log((N_NODES-1-noutedge-ninedge)/(noutedge+ninedge+1.0));
	  }else{
	    logratio = log((N_NODES-1-noutedge-ninedge)/(noutedge+ninedge+1.0));
	  }
      }
  }
  MHp->logratio += logratio;
}

/*********************
 void MH_ReallocateWithReplacement
*********************/
MH_P_FN(MH_ReallocateWithReplacement){  
  int i;
  Vertex root;
  Vertex* edges;
  int edgecount = 0;
  
  /* select a node at random */
  root = 1 + unif_rand() * N_NODES;

  edges = (Vertex *) Calloc(N_NODES+1, Vertex);
  for (i = 0; i <= N_NODES; i++)
    edges[i] = NO_EDGE;
  
  /* count current edges and mark them in an array */
  for (i = 1; i <= N_NODES; i++)
    {
      if (root == i) continue;
      if (IS_OUTEDGE(root, i) > 0)
	{
	  edges[i] = OLD_EDGE;
	  edgecount++;
	}
      if (!DIRECTED && (root > i) &&
	  (IS_OUTEDGE(i, root) > 0))
	{
	  edges[i] = OLD_EDGE;
	  edgecount++;
	}
    }
  
  /* select edgecount edges to create */
  for (i = 0; i < edgecount; i++)
    {
      Vertex newhead;
      /* get a new edge, neither the root nor something already chosen */
      while ((newhead = 1 + unif_rand() * N_NODES) == root ||
	     (edges[newhead] & NEW_EDGE))
	;
      
      /* if this edge already exists - (OLD_EDGE | NEW_EDGE) == CAN_IGNORE */
      edges[newhead] = edges[newhead] | NEW_EDGE;
    }
  
  /* index into Mtail/Mhead is  */
  edgecount = 0;
  
  /* add to toggle list:  anything that is non zero in edges array
     should be toggled, whether on or off. */
  for (i = 0; i <= N_NODES; i++)
    {
      if (edges[i] == NO_EDGE || edges[i] == CAN_IGNORE) continue;
      
      /* double to integer coercion */
      Mtail[edgecount] = root;
      Mhead[edgecount] = i;
      
      if (!DIRECTED && (Mtail[edgecount] > Mhead[edgecount]))
	{
	  Vertex temp;
	  temp = Mtail[edgecount];
	  Mtail[edgecount] = Mhead[edgecount];
	  Mhead[edgecount] = temp;
	}
      edgecount++;
    }
  Free(edges);
}

/*********************
 void MH_AllTogglesForOneNode
*********************/
MH_P_FN(MH_AllTogglesForOneNode){
  
  int i;
  int j;
  int root;
  
  root = 1 + unif_rand() * N_NODES;
  
  j = 0;
  for (i = 1; i <= N_NODES; i++)
    {
      /* probability here only do this with .8? */
      
      /* there is never an edge (root, root) */
      if (i == root)
	continue;
      
      /* double to integer coercion */
      Mtail[j] = root;
      Mhead[j] = i;
      
      if (!DIRECTED && (Mtail[j] > Mhead[j]))
	{
	  Vertex temp;
	  temp = Mtail[j];
	  Mtail[j] = Mhead[j];
	  Mhead[j] = temp;
	}
      j++;
    }
}


/*********************
 void MH_SwitchLabelTwoNodesToggles
*********************/
MH_P_FN(MH_SwitchLabelTwoNodesToggles){  
  int nedge1=0, nedge2=0, k, ntoggles;
  Vertex *edges1, *edges2;
  Vertex e, tail2, head2, tail1, head1;

  /* *** don't forget tail-> head now */
  
  /* select a node at random */
  edges1 = (Vertex *) Calloc(N_NODES+1, Vertex);
  edges2 = (Vertex *) Calloc(N_NODES+1, Vertex);
  
  while(nedge1==0){
    tail1 = 1 + unif_rand() * N_NODES;
    
    for(e = EdgetreeMinimum(nwp->outedges, tail1);
	(head1 = nwp->outedges[e].value) != 0; /* loop if */
	e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail1 */
      {
        edges1[nedge1] = head1;
	++nedge1;
      }
    for(e = EdgetreeMinimum(nwp->inedges, tail1);
	(head1 = nwp->inedges[e].value) != 0; /* loop if */
	e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail1 */
      {
        edges1[nedge1] = head1;
	++nedge1;
      }
  }
  
  while((tail2 = 1 + unif_rand() * N_NODES) == tail1);
  
  for(e = EdgetreeMinimum(nwp->outedges, tail2);
      (head2 = nwp->outedges[e].value) != 0; /* loop if */
      e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail2 */
    {
      edges2[nedge2] = head2;
      ++nedge2;
    }
  for(e = EdgetreeMinimum(nwp->inedges, tail2);
      (head2 = nwp->inedges[e].value) != 0; /* loop if */
      e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail2 */
    {
      edges2[nedge2] = head2;
      ++nedge2;
    }
  
  ntoggles = 0;
  for(k=0; k < nedge1; k++){
    if (tail1 > edges1[k])
      {
	Mtail[ntoggles] = edges1[k];
	Mhead[ntoggles] = tail1;
      }
    if (tail1 < edges1[k]){
      Mtail[ntoggles] = tail1;
      Mhead[ntoggles] = edges1[k];
    }
    if(tail1 != edges1[k]) ntoggles++;
  }
  
  for(k=0; k < nedge2; k++){
    if (tail1 > edges2[k])
      {
	Mtail[ntoggles] = edges2[k];
	Mhead[ntoggles] = tail1;
      }
    if (tail1 < edges2[k]){
      Mtail[ntoggles] = tail1;
      Mhead[ntoggles] = edges2[k];
    }
    if(tail1 != edges2[k]) ntoggles++;
  }
  
  for(k=0; k < nedge2; k++){
    if (tail2 > edges2[k])
      {
	Mtail[ntoggles] = edges2[k];
	Mhead[ntoggles] = tail2;
      }
    if (tail2 < edges2[k]){
      Mtail[ntoggles] = tail2;
      Mhead[ntoggles] = edges2[k];
    }
    if(tail2 != edges2[k]) ntoggles++;
  }
  
  for(k=0; k < nedge1; k++){
    if (tail2 > edges1[k])
      {
	Mtail[ntoggles] = edges1[k];
	Mhead[ntoggles] = tail2;
      }
    if (tail2 < edges1[k]){
      Mtail[ntoggles] = tail2;
      Mhead[ntoggles] = edges1[k];
    }
    if(tail2 != edges1[k]) ntoggles++;
  }
  Free(edges1);
  Free(edges2);
}


/*********************
 void MH_ConstrainedCondDegDist
*********************/
MH_P_FN(MH_ConstrainedCondDegDist){  
  int noutedge=0, ninedge=0, k, fvalid=0;
  int k0, j0, j1, k1;
  int j0h, j1h;
  Vertex *outedges, *inedges;
  Vertex e, alter, tail=0, head;

  /* *** don't forget tail-> head now */
  
  /* select a node at random */
  outedges = (Vertex *) Calloc(N_NODES+1, Vertex);
  inedges = (Vertex *) Calloc(N_NODES+1, Vertex);
  
  while(noutedge==0 && ninedge==0){
    tail = 1 + unif_rand() * N_NODES;
    
    for(e = EdgetreeMinimum(nwp->outedges, tail);
	(head = nwp->outedges[e].value) != 0; /* loop if */
	e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
      {
        outedges[noutedge] = head;
	++noutedge;
      }
    for(e = EdgetreeMinimum(nwp->inedges, tail);
	(head = nwp->inedges[e].value) != 0; /* loop if */
	e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail */
      {
        inedges[ninedge] = head;
	++ninedge;
      }
  }
  
  k0 = (int)(unif_rand() * (noutedge+ninedge)); 
  if (k0 < noutedge){
    head = outedges[k0]; 
  }else{
    head = inedges[k0-noutedge]; 
  }
  if ( (!DIRECTED && tail > head) ||
       (  DIRECTED  && k0 >= noutedge) )
    {
      Mtail[0] = head;
      Mhead[0] = tail;
    }else{
      Mtail[0] = tail;
      Mhead[0] = head;
    }
  
  if (dEdgeListSearch(Mtail[0], Mhead[0], MH_INPUTS)==0){
    Mtail[0] = Mhead[0] = 0;
    Mtail[1] = Mhead[1] = 0;
  }
  
  fvalid=0;
  k1=0;
  while(fvalid==0 && k1 < 100){
    while((alter = 1 + unif_rand() * N_NODES) == tail);
    if(alter != head){fvalid=1;}
    fvalid=1;
    if (k0 < noutedge || !(DIRECTED)){
      k=0;
      while(fvalid==1 && noutedge > 0 && k <= noutedge-1){
	if(alter == outedges[k]){fvalid=0;}else{++k;}
      }
    }
    if (k0 >= noutedge || !(DIRECTED)){
      k=0;
      while(fvalid==1 && ninedge > 0 && k <= ninedge-1){
	if(alter == inedges[k]){fvalid=0;}else{++k;}
      }
    }
    k1++;
  }
  
  if (k1 == 100){
    Mtail[0] = Mhead[0] = 0;
    Mtail[1] = Mhead[1] = 0;
  }
  
  if ( (!DIRECTED && alter > tail) ||
       (DIRECTED && k0 < noutedge) )
    {
      Mtail[1] = tail;
      Mhead[1] = alter;
    }else{
      Mtail[1] = alter;
      Mhead[1] = tail;
    }
  
  if (dEdgeListSearch(Mtail[1], Mhead[1], MH_INPUTS)==0){
    Mtail[0] = Mhead[0] = 0;
    Mtail[1] = Mhead[1] = 0;
  }
  
  Free(outedges);
  Free(inedges);
  
  /* Check undirected degrees */

  /* *** don't forget tail-> head now */

  if (!DIRECTED){
    k0=nwp->outdegree[tail]+ nwp->indegree[tail];
    j0h=nwp->outdegree[head]+ nwp->indegree[head];
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
      j0h=nwp->indegree[head];
      j1h=nwp->indegree[alter];
    }else{
      /* Check outdegrees */
      j0h=nwp->outdegree[head];
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
    Mtail[0] = Mhead[0] = 0;
    Mtail[1] = Mhead[1] = 0;
  }
}

void MH_ConstrainedNodePairedTiesToggles (MHProposal *MHp,
       	 Network *nwp) {  
  /* chooses a node and toggles all ties and
     and toggles an equal number of matching nonties
     for that node */
  int nedge=0,j,k;
  int fvalid = 1;
  Vertex e, tail, prop;

  /* *** don't forget tail-> head now */
  
  /* double to integer coercion */
  tail = 1 + unif_rand() * N_NODES; 
  
  for(e = EdgetreeMinimum(nwp->outedges, tail);
      (prop = nwp->outedges[e].value) != 0; /* loop if */
      e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
    {
      Mtail[nedge] = tail;
      Mhead[nedge] = prop;
      ++nedge;
    }
  for(e = EdgetreeMinimum(nwp->inedges, tail);
      (prop = nwp->inedges[e].value) != 0; /* loop if */
      e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail */
    {
      Mhead[nedge] = tail;
      Mtail[nedge] = prop;
      ++nedge;
    }
  
  if(nedge > N_NODES-nedge){
    Mtail[0] = Mhead[0] = 0;
    Mtail[1] = Mhead[1] = 0;
  }  
  j = 0;
  while (j <=nedge)
    {
      prop = 1 + unif_rand() * N_NODES; 
      k=0;
      fvalid=1;
      while(fvalid==1 && k<nedge+j){
	if(IS_OUTEDGE(MIN(prop,Mtail[k]),
			   MAX(prop,Mtail[k])) +
	   IS_OUTEDGE( MIN(prop,Mhead[k]),
			   MAX(prop,Mhead[k]))==0
	   ){++k;
	}else{
	  fvalid=0;}
      }
      if(prop>tail){
	Mtail[j+nedge] = tail;
	Mhead[j+nedge] = prop;
      }else{
	Mtail[j+nedge] = prop;
	Mhead[j+nedge] = tail;
      }
      ++j;
    }
  
  j = 2*nedge;
  if (!CheckConstrainedTogglesValid(MHp, nwp))
    {
      *Mtail = *Mhead = 0;
    }
}

/*********************
 void MH_ConstrainedReallocateWithReplacement
*********************/
void MH_ConstrainedReallocateWithReplacement (MHProposal *MHp,
       	 Network *nwp) {  
  int i;
  Vertex root;
  Vertex* edges;
  int edgecount = 0;
  
  /* select a node at random */
  root = 1 + unif_rand() * N_NODES;

  edges = (Vertex *) Calloc(N_NODES+1, Vertex);
  for (i = 0; i <= N_NODES; i++)
    edges[i] = NO_EDGE;
  
  /* count current edges and mark them in an array */
  for (i = 1; i <= N_NODES; i++)
    {
      if (root == i) continue;
      if (IS_OUTEDGE(root, i) > 0)
	{
	  edges[i] = OLD_EDGE;
	  edgecount++;
	}
      if (!DIRECTED && (root > i) &&
	  (IS_OUTEDGE(i, root) > 0))
	{
	  edges[i] = OLD_EDGE;
	  edgecount++;
	}
    }
  
  /* select edgecount edges to create */
  for (i = 0; i < edgecount; i++)
    {
      Vertex newhead;
      
      /* get a new edge, neither the root nor something already chosen */
      while ((newhead = 1 + unif_rand() * N_NODES) == root ||
	     (edges[newhead] & NEW_EDGE))
	;
      
      /* if this edge already exists - (OLD_EDGE | NEW_EDGE) == CAN_IGNORE */
      edges[newhead] = edges[newhead] | NEW_EDGE;
    }
  
  /* index into Mtail/Mhead is  */
  edgecount = 0;
  
  /* add to toggle list:  anything that is non zero in edges array
     should be toggled, whether on or off. */
  for (i = 0; i <= N_NODES; i++)
    {
      if (edges[i] == NO_EDGE || edges[i] == CAN_IGNORE) continue;
      
      /* double to integer coercion */
      Mtail[edgecount] = root;
      Mhead[edgecount] = i;
      
      if (!DIRECTED && (Mtail[edgecount] > Mhead[edgecount]))
	{
	  Vertex temp;
	  temp = Mtail[edgecount];
	  Mtail[edgecount] = Mhead[edgecount];
	  Mhead[edgecount] = temp;
	}
      edgecount++;
    }
  Free(edges);
}

/*********************
 void MH_ConstrainedAllTogglesForOneNode
*********************/
void MH_ConstrainedAllTogglesForOneNode (MHProposal *MHp,
					 Network *nwp) {
  int i;
  int j;
  int root;
  
  root = 1 + unif_rand() * N_NODES;
  
  j = 0;
  for (i = 1; i <= N_NODES; i++)
    {
      /* probability here only do this with .8? */
      
      /* there is never an edge (root, root) */
      if (i == root)
	continue;
      
      /* double to integer coercion */
      Mtail[j] = root;
      Mhead[j] = i;
      
      if (!DIRECTED && (Mtail[j] > Mhead[j]))
	{
	  Vertex temp;
	  temp = Mtail[j];
	  Mtail[j] = Mhead[j];
	  Mhead[j] = temp;
	}
      j++;
    }
}

/*********************
 void MH_ConstrainedTwoRandomToggles
*********************/
void MH_ConstrainedTwoRandomToggles (MHProposal *MHp,
				 Network *nwp) {  
  int i;
  
  for (i = 0; i < 2; i++)
    {
      /* double to integer coercion */
      Mtail[i] = 1 + unif_rand() * N_NODES; 
      while ((Mhead[i] = 1 + unif_rand() * N_NODES) == Mtail[i]);
      
      while(dEdgeListSearch(Mtail[i], Mhead[i], MH_INPUTS)==0){
	Mtail[i] = 1 + unif_rand() * N_NODES; 
	while ((Mhead[i] = 1 + unif_rand() * N_NODES) == Mtail[i]);
      }
      if (!DIRECTED && Mtail[i] > Mhead[i]) 
	{
	  Vertex temp;
	  temp = Mtail[i];
	  Mtail[i] = Mhead[i];
	  Mhead[i] = temp;
	}
    }
  
  if (!CheckConstrainedTogglesValid(MHp, nwp))
    {
      Mtail[0] = Mhead[0] = 0;
      Mtail[1] = Mhead[1] = 0;
    }  
}

/*********************
 void MH_ConstrainedCondDeg
*********************/
void MH_ConstrainedCondDeg (MHProposal *MHp,
					 Network *nwp) {  
  /* WARNING: THIS NEEDS TO BE FIXED */
  int nedge1=0, nedge2=0, k, toomany, fvalid=0;
  Vertex *edges1, *edges2;
  Vertex e, tail2=0, head2, tail1, head1;
  
  /* select a node at random */
  edges1 = (Vertex *) Calloc(N_NODES+1, Vertex);
  edges2 = (Vertex *) Calloc(N_NODES+1, Vertex);
  
  while(nedge1==0){
    tail1 = 1 + unif_rand() * N_NODES;
    
    for(e = EdgetreeMinimum(nwp->outedges, tail1);
	(head1 = nwp->outedges[e].value) != 0; /* loop if */
	e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail1 */
      {
        edges1[nedge1] = head1;
	++nedge1;
      }
    for(e = EdgetreeMinimum(nwp->inedges, tail1);
	(head1 = nwp->inedges[e].value) != 0; /* loop if */
	e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail1 */
      {
        edges1[nedge1] = head1;
	++nedge1;
      }
  }
  
  head1 = edges1[(int)(unif_rand() * nedge1)]; 
  if (tail1 > head1)
    {
      Mtail[0] = head1;
      Mhead[0] = tail1;
    }else{
      Mtail[0] = tail1;
      Mhead[0] = head1;
    }
   
  toomany = 0;
  while(nedge2==0 && toomany < 100){
    fvalid=0;
    while(fvalid==0){
      while((tail2 = 1 + unif_rand() * N_NODES) == tail1);
      k=0;
      fvalid=1;
      while(fvalid==1 && k < nedge1){
	if(tail2 == edges1[k]){fvalid=0;}else{++k;}
      }
    }

    for(e = EdgetreeMinimum(nwp->outedges, tail2);
	(head2 = nwp->outedges[e].value) != 0; /* loop if */
	e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail2 */
      {
        edges2[nedge2] = head2;
	++nedge2;
      }
    for(e = EdgetreeMinimum(nwp->inedges, tail2);
	(head2 = nwp->inedges[e].value) != 0; /* loop if */
	e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail2 */
      {
        edges2[nedge2] = head2;
	++nedge2;
      }
    ++toomany;
  }
  if (toomany==100){
    Mtail[0] = Mhead[0] = 0;
    Mtail[1] = Mhead[1] = 0;
  }
  toomany=0;
  fvalid=0;
  while(fvalid==0 && toomany < 10){
    while((head2 = edges2[(int)(unif_rand() * nedge2)]) == tail1);
    k=0;
    fvalid=1;
    while(fvalid==1 && k < nedge1){
      if(head2 == edges1[k]){fvalid=0;}else{++k;}
    }
    ++toomany;
  }
  if (!fvalid || toomany==10){
    Mtail[0] = Mhead[0] = 0;
    Mtail[1] = Mhead[1] = 0;
    Free(edges1);
    Free(edges2);
      }
  if (tail2 > head2)
    {
      Mtail[1] = head2;
      Mhead[1] = tail2;
    }else{
      Mtail[1] = tail2;
      Mhead[1] = head2;
    }
  Free(edges1);
  Free(edges2);
}

/*********************
 void MH_ConstrainedSwitchLabelTwoNodesToggles
*********************/
void MH_ConstrainedSwitchLabelTwoNodesToggles (MHProposal *MHp,
       	 Network *nwp)  {  
  int nedge1=0, nedge2=0, k, ntoggles;
  Vertex *edges1, *edges2;
  Vertex e, tail2, head2, tail1, head1;

  /* *** don't forget tail-> head now */
  
  /* select a node at random */

  edges1 = (Vertex *) Calloc(N_NODES+1, Vertex);
  edges2 = (Vertex *) Calloc(N_NODES+1, Vertex);

  while(nedge1==0){
    tail1 = 1 + unif_rand() * N_NODES;
    
    for(e = EdgetreeMinimum(nwp->outedges, tail1);
	(head1 = nwp->outedges[e].value) != 0; /* loop if */
	e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail1 */
      {
        edges1[nedge1] = head1;
	++nedge1;
      }
    for(e = EdgetreeMinimum(nwp->inedges, tail1);
	(head1 = nwp->inedges[e].value) != 0; /* loop if */
	e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail1 */
      {
        edges1[nedge1] = head1;
	++nedge1;
      }
  }
  
  while((tail2 = 1 + unif_rand() * N_NODES) == tail1);
  
  for(e = EdgetreeMinimum(nwp->outedges, tail2);
      (head2 = nwp->outedges[e].value) != 0; /* loop if */
      e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail2 */
    {
      edges2[nedge2] = head2;
      ++nedge2;
    }
  for(e = EdgetreeMinimum(nwp->inedges, tail2);
      (head2 = nwp->inedges[e].value) != 0; /* loop if */
      e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail2 */
    {
      edges2[nedge2] = head2;
      ++nedge2;
    }
  
  ntoggles = 0;
  for(k=0; k < nedge1; k++){
    if (tail1 > edges1[k])
      {
	Mtail[ntoggles] = edges1[k];
	Mhead[ntoggles] = tail1;
      }
    if (tail1 < edges1[k]){
      Mtail[ntoggles] = tail1;
      Mhead[ntoggles] = edges1[k];
    }
    if(tail1 != edges1[k]) ntoggles++;
  }
  
  for(k=0; k < nedge2; k++){
    if (tail1 > edges2[k])
      {
	Mtail[ntoggles] = edges2[k];
	Mhead[ntoggles] = tail1;
      }
    if (tail1 < edges2[k]){
      Mtail[ntoggles] = tail1;
      Mhead[ntoggles] = edges2[k];
    }
    if(tail1 != edges2[k]) ntoggles++;
  }
  
  for(k=0; k < nedge2; k++){
    if (tail2 > edges2[k])
      {
	Mtail[ntoggles] = edges2[k];
	Mhead[ntoggles] = tail2;
      }
    if (tail2 < edges2[k]){
      Mtail[ntoggles] = tail2;
      Mhead[ntoggles] = edges2[k];
    }
    if(tail2 != edges2[k]) ntoggles++;
  }
  
  for(k=0; k < nedge1; k++){
    if (tail2 > edges1[k])
      {
	Mtail[ntoggles] = edges1[k];
	Mhead[ntoggles] = tail2;
      }
    if (tail2 < edges1[k]){
      Mtail[ntoggles] = tail2;
      Mhead[ntoggles] = edges1[k];
    }
    if(tail2 != edges1[k]) ntoggles++;
  }
  Free(edges1);
  Free(edges2);
}

/*********************
 void MH_ConstantEdgesToggles
*********************/
MH_P_FN(MH_ConstantEdgesToggles){  
  int noutedge=0, ninedge=0, k, fvalid=0;
  int k0, k1;
  Vertex e, alter, tail, head, head1;

  /* *** don't forget tail-> head now */
  
  while(noutedge+ninedge==0){
    /* select a node at random */
    tail = 1 + unif_rand() * N_NODES;
    ninedge  = nwp->indegree[tail];
    noutedge = nwp->outdegree[tail];
  }
  
  k0 = (int)(unif_rand() * (noutedge+ninedge)); 
  if (k0 < noutedge){
    k=0;
    for(e = EdgetreeMinimum(nwp->outedges, tail);
	((head = nwp->outedges[e].value) != 0 && k<k0);
	e = EdgetreeSuccessor(nwp->outedges, e)){++k;}
  }else{
    k=0;
    for(e = EdgetreeMinimum(nwp->inedges, tail);
	((head = nwp->inedges[e].value) != 0 && k<(k0-noutedge));
	e = EdgetreeSuccessor(nwp->inedges, e)){++k;}
  }
  
  if ( (!DIRECTED && tail > head) ||
       (DIRECTED && k0 >= noutedge) )
    {
      Mtail[0] = head;
      Mhead[0] = tail;
    }else{
      Mtail[0] = tail;
      Mhead[0] = head;
    }
  
  k1=0;
  fvalid=0;
  while(fvalid==0 && k1 < 100){
    while((alter = 1 + unif_rand() * N_NODES) == tail);
    fvalid=1;
    if(alter == head){fvalid=0;}
    if (k0 < noutedge || !(DIRECTED)){
      for(e = EdgetreeMinimum(nwp->outedges, tail);
	  (fvalid==1 && ((head1 = nwp->outedges[e].value) != 0));
	  e = EdgetreeSuccessor(nwp->outedges, e)){
	if(alter==head1){fvalid=0;}}
    }
    if (k0 >= noutedge || !(DIRECTED)){
      for(e = EdgetreeMinimum(nwp->inedges, tail);
	  (fvalid==1 && ((head1 = nwp->inedges[e].value) != 0));
	  e = EdgetreeSuccessor(nwp->inedges, e)){
	if(alter==head1){fvalid=0;}}
    }
    k1++;
  }
  if (k1 == 100){
    Mtail[0] = Mhead[0] = 0;
    Mtail[1] = Mhead[1] = 0;
  }
  
  if ( (!DIRECTED && alter > tail) ||
       (DIRECTED && k0 < noutedge) )
    {
      Mtail[1] = tail;
      Mhead[1] = alter;
    }else{
      Mtail[1] = alter;
      Mhead[1] = tail;
    }
  
  if (!fvalid){
    Mtail[0] = Mhead[0] = 0;
    Mtail[1] = Mhead[1] = 0;
  }else{  
  }
}

/*********************
 void MH_CondDegSwitchToggles
*********************/
MH_P_FN(MH_CondDegSwitchToggles){  
  int noutedge, ninedge, i;
  int k, k0, toomany;
  Vertex e, tail, head;

  /* *** don't forget tail-> head now */
  
  /* select a node at random */
  for (i = 0; i < 2; i++){
    toomany=0;
    noutedge=0;
    ninedge=0;
    while(noutedge==0 && ninedge==0 && toomany < 100){
      tail = 1 + unif_rand() * N_NODES;
      ninedge=0;
      noutedge=0;
      while(noutedge+ninedge==0){
	/* select a node at random */
	tail = 1 + unif_rand() * N_NODES;
	ninedge = nwp->indegree[tail];
	noutedge = nwp->outdegree[tail];
      }
      ++toomany;
    }
    
    if (toomany == 100){
      Mtail[0] = Mhead[0] = 0;
      Mtail[1] = Mhead[1] = 0;
    }
    
    k0 = (int)(unif_rand() * (noutedge+ninedge)); 
    if (k0 < noutedge){
      k=0;
      for(e = EdgetreeMinimum(nwp->outedges, tail);
	  ((head = nwp->outedges[e].value) != 0 && k<k0);
	  e = EdgetreeSuccessor(nwp->outedges, e)){++k;}
    }else{
      k=0;
      for(e = EdgetreeMinimum(nwp->inedges, tail);
	  ((head = nwp->inedges[e].value) != 0 && k<(k0-noutedge));
	  e = EdgetreeSuccessor(nwp->inedges, e)){++k;}
    }
    if ( (!DIRECTED && tail > head) ||
	 (DIRECTED && k0 >= noutedge) )
      {
	Mtail[i] = head;
	Mhead[i] = tail;
      }else{
	Mtail[i] = tail;
	Mhead[i] = head;
      }
  }
  
  if (IS_OUTEDGE( Mtail[0],Mhead[1]) ||
      IS_OUTEDGE( Mtail[1],Mhead[0]) ){
    Mtail[0] = Mhead[0] = 0;
    Mtail[1] = Mhead[1] = 0;
  }
  
  if ( (!DIRECTED && Mtail[0] > Mhead[1]) )
    {
      Mtail[2] = Mhead[1];
      Mhead[2] = Mtail[0];
    }else{
      Mtail[2] = Mtail[0];
      Mhead[2] = Mhead[1];
    }
  
  if ( (!DIRECTED && Mtail[1] > Mhead[0]) )
    {
      Mtail[3] = Mhead[0];
      Mhead[3] = Mtail[1];
    }else{
      Mtail[3] = Mtail[1];
      Mhead[3] = Mhead[0];
    }
}


