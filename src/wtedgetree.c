/*  File src/wtedgetree.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2013 Statnet Commons
 */
#include "wtedgetree.h"

/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */


/*******************
 WtNetwork WtNetworkInitialize

 Initialize, construct binary tree version of network with weights.  Note
 that the 0th WtTreeNode in the array is unused and should 
 have all its values set to zero
*******************/
/* *** don't forget, tail -> head */

WtNetwork WtNetworkInitialize(Vertex *tails, Vertex *heads, double *weights,
			      Edge nedges, Vertex nnodes, int directed_flag, Vertex bipartite,
			      int lasttoggle_flag, int time, int *lasttoggle, unsigned int n_aux) {
  WtNetwork nw;

  nw.last_inedge = nw.last_outedge = (Edge)nnodes;
  /* Calloc will zero the allocated memory for us, probably a lot
     faster. */
  nw.outdegree = (Vertex *) calloc((nnodes+1),sizeof(Vertex));
  nw.indegree  = (Vertex *) calloc((nnodes+1),sizeof(Vertex));
  nw.maxedges = MAX(nedges,1)+nnodes+2; /* Maybe larger than needed? */
  nw.inedges = (WtTreeNode *) calloc(nw.maxedges,sizeof(WtTreeNode));
  nw.outedges = (WtTreeNode *) calloc(nw.maxedges,sizeof(WtTreeNode));

  GetRNGstate();  /* R function enabling uniform RNG */

  if(lasttoggle_flag){
    nw.duration_info.time=time;
    if(lasttoggle){
      nw.duration_info.lasttoggle = (int *) calloc(DYADCOUNT(nnodes, bipartite, directed_flag), sizeof(int));
      memcpy(nw.duration_info.lasttoggle, lasttoggle, DYADCOUNT(nnodes, bipartite, directed_flag) * sizeof(int));
    } else nw.duration_info.lasttoggle = NULL;
  }
  else nw.duration_info.lasttoggle = NULL;

  /*Configure a Network*/
  nw.nnodes = nnodes;
  nw.nedges = 0; /* Edges will be added one by one */
  nw.directed_flag=directed_flag;
  nw.bipartite=bipartite;

  WtShuffleEdges(tails,heads,weights,nedges); /* shuffle to avoid worst-case performance */

  for(Edge i = 0; i < nedges; i++) {
    Vertex tail=tails[i], head=heads[i];
    double w=weights[i];
    if (!directed_flag && tail > head) 
      WtAddEdgeToTrees(head,tail,w,&nw); /* Undir edges always have tail < head */ 
    else 
      WtAddEdgeToTrees(tail,head,w,&nw);
  }

  /* Allocate pointers to auxiliary storage */
  nw.n_aux = n_aux;
  if(n_aux){
    nw.aux_storage = (void **)malloc(sizeof(void *)*n_aux);
    for(unsigned int i = 0; i<n_aux; i++) nw.aux_storage[i] = NULL;
  }else nw.aux_storage = NULL;
  
  PutRNGstate();  
  return nw;
}


/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */

/*Takes vectors of doubles for edges; used only when constructing from inputparams. */
WtNetwork WtNetworkInitializeD(double *tails, double *heads, double *weights, Edge nedges,
           Vertex nnodes, int directed_flag, Vertex bipartite,
           int lasttoggle_flag, int time, int *lasttoggle, unsigned int n_aux) {

  /* *** don't forget, tail -> head */

  Vertex *itails=(Vertex*)malloc(sizeof(Vertex)*nedges);
  Vertex *iheads=(Vertex*)malloc(sizeof(Vertex)*nedges);
  
  for(Edge i=0; i<nedges; i++){
    itails[i]=tails[i];
    iheads[i]=heads[i];
  }

  WtNetwork nw=WtNetworkInitialize(itails,iheads,weights,nedges,nnodes,directed_flag,bipartite,lasttoggle_flag, time, lasttoggle, n_aux);

  free(itails);
  free(iheads);
  return nw;
}

/*******************
 void NetworkDestroy
*******************/
void WtNetworkDestroy(WtNetwork *nwp) {
  free (nwp->indegree);
  free (nwp->outdegree);
  free (nwp->inedges);
  free (nwp->outedges);
  if(nwp->duration_info.lasttoggle){
    free (nwp->duration_info.lasttoggle);
    nwp->duration_info.lasttoggle=NULL;
  }
  if(nwp->aux_storage){
    for(unsigned int i=0; i<nwp->n_aux; i++)
      if(nwp->aux_storage[i]){
	  free(nwp->aux_storage[i]);
	  nwp->aux_storage[i] = NULL;
      }
    free(nwp->aux_storage);
    nwp->aux_storage = NULL;
  }
}

/******************
 Network WtNetworkCopy
*****************/
WtNetwork *WtNetworkCopy(WtNetwork *dest, WtNetwork *src){
  Vertex nnodes = dest->nnodes = src->nnodes;
  dest->last_inedge = src->last_inedge;
  dest->last_outedge = src->last_outedge;

  dest->outdegree = (Vertex *) malloc((nnodes+1)*sizeof(Vertex));
  memcpy(dest->outdegree, src->outdegree, (nnodes+1)*sizeof(Vertex));
  dest->indegree = (Vertex *) malloc((nnodes+1)*sizeof(Vertex));
  memcpy(dest->indegree, src->indegree, (nnodes+1)*sizeof(Vertex));

  Vertex maxedges = dest->maxedges = src->maxedges;

  dest->inedges = (WtTreeNode *) malloc(maxedges*sizeof(WtTreeNode));
  memcpy(dest->inedges, src->inedges, maxedges*sizeof(WtTreeNode));
  dest->outedges = (WtTreeNode *) malloc(maxedges*sizeof(WtTreeNode));
  memcpy(dest->outedges, src->outedges, maxedges*sizeof(WtTreeNode));

  int directed_flag = dest->directed_flag = src->directed_flag;
  Vertex bipartite = dest->bipartite = src->bipartite;

  if(src->duration_info.lasttoggle){
    dest->duration_info.time=src->duration_info.time;
    dest->duration_info.lasttoggle = (int *) calloc(DYADCOUNT(nnodes, bipartite, directed_flag), sizeof(int));
    memcpy(dest->duration_info.lasttoggle, src->duration_info.lasttoggle,DYADCOUNT(nnodes, bipartite, directed_flag) * sizeof(int));
  }
  else dest->duration_info.lasttoggle = NULL;

  dest->nedges = src->nedges;

  /* Allocate pointers to auxiliary storage and set them to NULL. Change stats should know to reinitialize. */
  if(src->n_aux){
    dest->n_aux = src->n_aux;
    dest->aux_storage = (void **)malloc(sizeof(void *)*dest->n_aux);
    for(unsigned int i = 0; i<dest->n_aux; i++) dest->aux_storage[i] = NULL;
  }
  
  return dest;
}

/*****************
 void Wtprintedge

 Diagnostic routine that prints out the contents
 of the specified WtTreeNode (used for debugging).  
*****************/
void Wtprintedge(Edge e, WtTreeNode *edges){
  Rprintf("Edge structure [%d]:\n",e);
  Rprintf("\t.value=%d\n",edges[e].value);
  Rprintf("\t.parent=%d\n",edges[e].parent);
  Rprintf("\t.left=%d\n",edges[e].left);
  Rprintf("\t.right=%d\n",edges[e].right);
  Rprintf("\t.weight=%d\n",edges[e].weight);
}

/*****************
 void WtNetworkEdgeList

 Print an edgelist for a network
*****************/
void WtNetworkEdgeList(WtNetwork *nwp) {
  Vertex i;
  for(i=1; i<=nwp->nnodes; i++) {
    Rprintf("Node %d:\n  ", i);
    WtInOrderTreeWalk(nwp->outedges, i);
    Rprintf("\n");
  }
}

/*****************
 void WtInOrderTreeWalk

 Diagnostic routine that prints the nodes in the tree rooted
 at edges[x], in increasing order of their values.
*****************/
void WtInOrderTreeWalk(WtTreeNode *edges, Edge x) {
  if (x != 0) {
    WtInOrderTreeWalk(edges, (edges+x)->left);
    /*    printedge(x, edges); */
    Rprintf(" %d:%f ",(edges+x)->value, (edges+x)->weight); 
    WtInOrderTreeWalk(edges, (edges+x)->right);
  }
}

/*****************
 Edge DesignMissing (see EdgetreeSearch)
*****************/
Edge WtDesignMissing (Vertex a, Vertex b, WtNetwork *mnwp) {
  int miss;
  miss = WtEdgetreeSearch(a,b,mnwp->outedges);
  if(mnwp->directed_flag){
    miss += WtEdgetreeSearch(a,b,mnwp->inedges);
  }
  return(miss);
}


/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */


/*****************
  int WtFindithEdge

  Find the ith edge in the WtNetwork *nwp and update the values of
  tail, head, and weight appropriately. If the value passed to tail,
  head, or weight is NULL, it is not updated, so it is possible to
  only obtain what is needed. Return 1 if successful, 0 otherwise.
  Note that i is numbered from 1, not 0.  Thus, the maximum possible
  value of i is nwp->nedges.
******************/

/* *** don't forget tail->head, so this function now accepts tail before head */

int WtFindithEdge (Vertex *tail, Vertex *head, double *weight, Edge i, WtNetwork *nwp) {
  Vertex taili=1;
  Edge e;

  /* TODO: This could be speeded up by a factor of 3 or more by starting
     the search from the tail n rather than tail 1 if i > ndyads/2. */

  if (i > nwp->nedges || i<=0)
    return 0;
  while (i > nwp->outdegree[taili]) {
    i -= nwp->outdegree[taili];
    taili++;
  }

  /* TODO: This could be speeded up by a factor of 3 or more by starting
     the search from the tree maximum rather than minimum (left over) i > outdegree[taili]. */

  e=WtEdgetreeMinimum(nwp->outedges,taili);
  while (i-- > 1) {
    e=WtEdgetreeSuccessor(nwp->outedges, e);
  }
  if(tail) *tail = taili;
  if(head) *head = nwp->outedges[e].value;
  if(weight) *weight = nwp->outedges[e].weight;
  return 1;
}

/*****************
  int GetRandEdge

  Select an edge in the Network *nwp at random and update the values
  of tail and head appropriately. Return 1 if successful, 0 otherwise.
******************/

/* *** don't forget tail->head, so this function now accepts tail before head */

int WtGetRandEdge(Vertex *tail, Vertex *head, double *weight, WtNetwork *nwp) {
  if(nwp->nedges==0) return(0);
  // FIXME: The constant maxEattempts needs to be tuned.
  const unsigned int maxEattempts=10;
  unsigned int Eattempts = (nwp->maxedges-1)/nwp->nedges;
  Edge rane;
  
  if(Eattempts>maxEattempts){
    // If the outedges is too sparse, revert to the old algorithm.
    rane=1 + unif_rand() * nwp->nedges;
    WtFindithEdge(tail, head, weight, rane, nwp);
  }else{
    // Otherwise, find a TreeNode which has a head.
    do{
      // Note that the outedges array has maxedges elements, but the
      // 0th one is always blank, and those with index >
      // nwp->last_outedge are blank as well, so we need to generate
      // an index from 1 through nwp->last_outedge (inclusive).
      rane = 1 + unif_rand() * nwp->last_outedge;
    }while((*head=nwp->outedges[rane].value)==0); // Form the head, while we are at it.

    // Get the weight.
    if(weight)
      *weight=nwp->outedges[rane].weight;
    
    // Ascend the edgetree as long as we can.
    // Note that it will stop as soon as rane no longer has a parent,
    // _before_ overwriting it.
    while(nwp->outedges[rane].parent) rane=nwp->outedges[rane].parent;
    // Then, the position of the root is the tail of edge.
    *tail=rane;
  }
  return 1;
}

/*****************
  int WtFindithNonedge

  Find the ith nonedge in the WtNetwork *nwp and
  update the values of tail and head appropriately.  Return
  1 if successful, 0 otherwise.  
  Note that i is numbered from 1, not 0.  Thus, the maximum possible
  value of i is (ndyads - nwp->nedges).
******************/

/* This function is not yet written.  It's not clear whether it'll
   be needed. */      
  /* *** but if it is needed, don't forget,  tail -> head */

int WtFindithNonedge (Vertex *tail, Vertex *head, Dyad i, WtNetwork *nwp) {
  Vertex taili=1;
  Edge e;
  Dyad ndyads = DYADCOUNT(nwp->nnodes, nwp->bipartite, nwp->directed_flag);
  Vertex nheads = nwp->bipartite ? nwp->nnodes-nwp->bipartite : nwp->nnodes-1;
  
  // If the index is too high or too low, exit immediately.
  if (i > ndyads - nwp->nedges || i<=0)
    return 0;

  /* TODO: This could be speeded up by a factor of 3 or more by starting
     the search from the tail n rather than tail 1 if i > ndyads/2. */


  while (i > nheads - nwp->outdegree[taili]) {   // nheads - nwp->oudegree[taili] is the number of nonoutties of taili.
    i -= nheads - nwp->outdegree[taili];
    taili++;
  }

  // Now, our tail is taili.

  /* TODO: This could be speeded up by a factor of 3 or more by starting
     the search from the tree maximum rather than minimum (left over) i > outdegree[taili]. */

  Vertex lhead = 0;
  e = WtEdgetreeMinimum(nwp->outedges,taili);
  Vertex rhead = nwp->outedges[e].value;
  // Note that rhead-lhead-1 is the number of nonties between two successive ties.
  while (i > rhead-lhead-1) {
    i -= rhead-lhead-1;
    lhead = rhead;
    e = WtEdgetreeSuccessor(nwp->outedges, e);
    // If rhead was the highest-indexed head, then e is now 0.
    if(e) rhead = nwp->outedges[e].value;
    else break; // Note that we don't actually need rhead in the final step.
  }

  // Now, the head we are looking for is (left over) i after lhead.

  *tail = taili;
  *head = lhead + i;
  return 1;
}

/*****************
  int WtGetRandNonedge

  Select an non-edge in the WtNetwork *nwp at random and update the values
  of tail and head appropriately. Return 1 if successful, 0 otherwise.
******************/

/* *** don't forget tail->head, so this function now accepts tail before head */

int WtGetRandNonedge(Vertex *tail, Vertex *head, WtNetwork *nwp) {
  Dyad ndyads = DYADCOUNT(nwp->nnodes, nwp->bipartite, nwp->directed_flag);
  if(ndyads-nwp->nedges==0) return(0);

  /* There are two ways to get a random nonedge: 1) keep trying dyads
     at random until you find one that's not an edge or 2) generate i
     at random and find ith nonedge. Method 1 works better in sparse
     networks, while Method 2, which runs in deterministic time, works
     better in dense networks.

     The expected number of attempts for Method 1 is 1/(1-e/d) =
     d/(d-e), where e is the number of edges and d is the number of
     dyads.
  */

  // FIXME: The constant maxEattempts needs to be tuned.
  const unsigned int maxEattempts=10;
  unsigned int Eattempts = ndyads/(ndyads-nwp->nedges);
  
  if(Eattempts>maxEattempts){
    // If the network is too dense, use the deterministic-time method:
    Dyad rane=1 + unif_rand() * (ndyads-nwp->nedges);
    WtFindithNonedge(tail, head, rane, nwp);
  }else{
    do{
      GetRandDyad(tail, head, nwp);
    }while(WtEdgetreeSearch(*tail, *head, nwp->outedges));
  }
  return 1;
}

/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */


/****************
 Edge EdgeTree2EdgeList

 Write the edgelist of a network into tail and head arrays.
 Returns the number of edges in the network.
****************/
Edge WtEdgeTree2EdgeList(Vertex *tails, Vertex *heads, double *weights, WtNetwork *nwp, Edge nmax){
  Edge nextedge=0;

  /* *** don't forget,  tail -> head */
  if (nwp->directed_flag) {
    for (Vertex v=1; v<=nwp->nnodes; v++){
      for(Vertex e = WtEdgetreeMinimum(nwp->outedges, v);
      nwp->outedges[e].value != 0 && nextedge < nmax;
      e = WtEdgetreeSuccessor(nwp->outedges, e)){
        tails[nextedge] = v;
        heads[nextedge] = nwp->outedges[e].value;
	if(weights) weights[nextedge] = nwp->outedges[e].weight;
        nextedge++;
      }
    }
  }else{
    for (Vertex v=1; v<=nwp->nnodes; v++){
      for(Vertex e = WtEdgetreeMinimum(nwp->outedges, v);
      nwp->outedges[e].value != 0 && nextedge < nmax;
      e = WtEdgetreeSuccessor(nwp->outedges, e)){
        Vertex k = nwp->outedges[e].value;
        if(v < k){
          tails[nextedge] = k;
          heads[nextedge] = v;
	  if(weights) weights[nextedge] = nwp->outedges[e].weight;
          nextedge++;
        }else{
          tails[nextedge] = v;
          heads[nextedge] = k;
	  if(weights) weights[nextedge] = nwp->outedges[e].weight;
          nextedge++;
        }
      }
    }
  }
  return nextedge;
}


/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */


void WtShuffleEdges(Vertex *tails, Vertex *heads, double *weights, Edge nedges){
  /* *** don't forget,  tail -> head */
  for(Edge i = nedges; i > 0; i--) {
    Edge j = i * unif_rand();
    Vertex tail = tails[j];
    Vertex head = heads[j];
    double w = weights[j];
    tails[j] = tails[i-1];
    heads[j] = heads[i-1];
    weights[j] = weights[i-1];
    tails[i-1] = tail;
    heads[i-1] = head;
    weights[i-1] = w;
  }
}

#ifndef INLINE_EDGETREE
#include "wtedgetree_inline.inc"
#endif
