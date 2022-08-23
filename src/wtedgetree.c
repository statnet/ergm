/*  File src/wtedgetree.c in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2022 Statnet Commons
 */
#include "ergm_wtedgetree.h"
#include "ergm_Rutil.h"

#include "wtedgetree_inline.do_not_include_directly.h"

/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */


/*******************
 WtNetwork WtNetworkInitialize

 Initialize, construct binary tree version of network with weights.  Note
 that the 0th WtTreeNode in the array is unused and should 
 have all its values set to zero

Note: passing nedges > 0 and tails == heads == NULL is OK: it creates an empty network with at least nedges of space preallocated.
*******************/
/* *** don't forget, tail -> head */

WtNetwork *WtNetworkInitialize(Vertex *tails, Vertex *heads, double *weights,
			       Edge nedges, Vertex nnodes, int directed_flag, Vertex bipartite,
			       int lasttoggle_flag, int time, int *lasttoggle) {
  WtNetwork *nwp = Calloc(1, WtNetwork);

  nwp->eattrname = NULL;

  nwp->last_inedge = nwp->last_outedge = (Edge)nnodes;
  /* Calloc will zero the allocated memory for us, probably a lot
     faster. */
  nwp->outdegree = (Vertex *) Calloc((nnodes+1), Vertex);
  nwp->indegree  = (Vertex *) Calloc((nnodes+1), Vertex);
  nwp->maxedges = MAX(nedges,1)+nnodes+2; /* Maybe larger than needed? */
  nwp->inedges = (WtTreeNode *) Calloc(nwp->maxedges, WtTreeNode);
  nwp->outedges = (WtTreeNode *) Calloc(nwp->maxedges, WtTreeNode);

  if(lasttoggle_flag) error("The lasttoggle API has been removed from ergm.");

  /*Configure a Network*/
  nwp->nnodes = nnodes;
  EDGECOUNT(nwp) = 0; /* Edges will be added one by one */
  nwp->directed_flag=directed_flag;
  nwp->bipartite=bipartite;

  if(tails==NULL && heads==NULL && weights==NULL) return nwp;

  WtDetShuffleEdges(tails,heads,weights,nedges); /* shuffle to avoid worst-case performance */

  for(Edge i = 0; i < nedges; i++) {
    Vertex tail=tails[i], head=heads[i];
    double w=weights[i];
    if(w==0) continue;
    if (!directed_flag && tail > head) 
      WtAddEdgeToTrees(head,tail,w,nwp); /* Undir edges always have tail < head */ 
    else 
      WtAddEdgeToTrees(tail,head,w,nwp);
  }

  WtDetUnShuffleEdges(tails,heads,weights,nedges); /* Unshuffle edges */

  return nwp;
}


/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */

/*Takes vectors of doubles for edges; used only when constructing from inputparams. */
WtNetwork *WtNetworkInitializeD(double *tails, double *heads, double *weights, Edge nedges,
           Vertex nnodes, int directed_flag, Vertex bipartite,
           int lasttoggle_flag, int time, int *lasttoggle) {

  /* *** don't forget, tail -> head */

  Vertex *itails=(Vertex*)Calloc(nedges, Vertex);
  Vertex *iheads=(Vertex*)Calloc(nedges, Vertex);
  
  for(Edge i=0; i<nedges; i++){
    itails[i]=tails[i];
    iheads[i]=heads[i];
  }

  WtNetwork *nwp=WtNetworkInitialize(itails,iheads,weights,nedges,nnodes,directed_flag,bipartite,lasttoggle_flag, time, lasttoggle);

  Free(itails);
  Free(iheads);
  return nwp;
}

/*******************
 void NetworkDestroy
*******************/
void WtNetworkDestroy(WtNetwork *nwp) {
  Free(nwp->on_edge_change);
  Free(nwp->on_edge_change_payload);
  Free(nwp->indegree);
  Free(nwp->outdegree);
  Free(nwp->inedges);
  Free(nwp->outedges);
  Free(nwp);
}

/******************
 Network WtNetworkCopy
*****************/
WtNetwork *WtNetworkCopy(WtNetwork *src){
  WtNetwork *dest = Calloc(1, WtNetwork);

  Vertex nnodes = dest->nnodes = src->nnodes;
  dest->last_inedge = src->last_inedge;
  dest->last_outedge = src->last_outedge;

  dest->outdegree = (Vertex *) Calloc((nnodes+1), Vertex);
  memcpy(dest->outdegree, src->outdegree, (nnodes+1)*sizeof(Vertex));
  dest->indegree = (Vertex *) Calloc((nnodes+1), Vertex);
  memcpy(dest->indegree, src->indegree, (nnodes+1)*sizeof(Vertex));

  Vertex maxedges = dest->maxedges = src->maxedges;

  dest->inedges = (WtTreeNode *) Calloc(maxedges, WtTreeNode);
  memcpy(dest->inedges, src->inedges, maxedges*sizeof(WtTreeNode));
  dest->outedges = (WtTreeNode *) Calloc(maxedges, WtTreeNode);
  memcpy(dest->outedges, src->outedges, maxedges*sizeof(WtTreeNode));

  dest->directed_flag = src->directed_flag;
  dest->bipartite = src->bipartite;

  EDGECOUNT(dest) = EDGECOUNT(src);

  return dest;
}

/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */



/*****************
 int WtToggleEdge

 Toggle an edge:  Set it to the opposite of its current
 value.  Return 1 if edge added, 0 if deleted.
*****************/

/* *** don't forget tail->head, so this function now accepts tail before head */

int WtToggleEdge (Vertex tail, Vertex head, double weight, WtNetwork *nwp) 
{
  /* don't forget tails < heads now for undirected networks */
  ENSURE_TH_ORDER;
  if(WtDeleteEdgeFromTrees(tail,head,nwp))
    return 0;
  else{
    WtAddEdgeToTrees(tail,head,weight,nwp);
    return 1;
  }
}

/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */

/*****************
 Edge AddEdgeToTrees

 Add an edge from tail to head; note that the function assumes that it
 is legal and therefore does not have a return value.  Since each
 "edge" should be added to both the list of outedges and the list of
 inedges, this actually involves two calls to AddHalfedgeToTree (hence
 "Trees" instead of "Tree" in the name of this function).
*****************/

/* *** don't forget tail->head, so this function now accepts tail before head */

void WtAddEdgeToTrees(Vertex tail, Vertex head, double weight, WtNetwork *nwp){
#ifdef DEBUG
  if(WtEdgetreeSearch(tail, head, nwp->outedges)||WtEdgetreeSearch(head, tail, nwp->inedges)) error("WtAddEdgeToTrees() called for an extant edge. Note that this produces an error only if compiling with DEBUG macro set and silently produces undefined behavior otherwise.");
#endif // DEBUG
  for(unsigned int i = 0; i < nwp->n_on_edge_change; i++) nwp->on_edge_change[i](tail, head, weight, nwp->on_edge_change_payload[i], nwp, 0);

  WtAddHalfedgeToTree(tail, head, weight, nwp->outedges, &(nwp->last_outedge));
  WtAddHalfedgeToTree(head, tail, weight, nwp->inedges, &(nwp->last_inedge));
  ++nwp->outdegree[tail];
  ++nwp->indegree[head];
  ++EDGECOUNT(nwp);
  WtCheckEdgetreeFull(nwp);
}

/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */

/*****************
 int WtDeleteEdgeFromTrees

 Find and delete the edge from tail to head.  
 Return 1 if successful, 0 otherwise.  As with AddEdgeToTrees, this must
 be done once for outedges and once for inedges.
*****************/

/* *** don't forget tail->head, so this function now accepts tail before head */

int WtDeleteEdgeFromTrees(Vertex tail, Vertex head, WtNetwork *nwp){
  Edge zth, zht;
  if((zth=WtEdgetreeSearch(tail, head, nwp->outedges))&&(zht=WtEdgetreeSearch(head, tail, nwp->inedges))){
    if(nwp->n_on_edge_change){
      double w = nwp->outedges[zth].weight;
      for(unsigned int i = 0; i < nwp->n_on_edge_change; i++) nwp->on_edge_change[i](tail, head, 0, nwp->on_edge_change_payload[i], nwp, w);
    }
    WtDeleteHalfedgeFromTreeAt(tail, head, nwp->outedges,&(nwp->last_outedge), zth);
    WtDeleteHalfedgeFromTreeAt(head, tail, nwp->inedges, &(nwp->last_inedge), zht);
    --nwp->outdegree[tail];
    --nwp->indegree[head];
    --EDGECOUNT(nwp);
    return 1;
  }
  return 0;
}

/*****************
 void AddOnWtNetworkEdgeChange

 Insert a specified toggle callback at the specified position.
*****************/
void AddOnWtNetworkEdgeChange(WtNetwork *nwp, OnWtNetworkEdgeChange callback, void *payload, unsigned int pos){
  if(nwp->n_on_edge_change+1 > nwp->max_on_edge_change){
    nwp->max_on_edge_change = MAX(nwp->max_on_edge_change,1)*2;
    nwp->on_edge_change = Realloc(nwp->on_edge_change, nwp->max_on_edge_change, OnWtNetworkEdgeChange);
    nwp->on_edge_change_payload = Realloc(nwp->on_edge_change_payload, nwp->max_on_edge_change, void*);
  }

  pos = MIN(pos, nwp->n_on_edge_change); // Last position.
  // Move everything down the list.
  memmove(nwp->on_edge_change+pos+1, nwp->on_edge_change+pos, (nwp->n_on_edge_change-pos)*sizeof(OnWtNetworkEdgeChange));
  memmove(nwp->on_edge_change_payload+pos+1, nwp->on_edge_change_payload+pos, (nwp->n_on_edge_change-pos)*sizeof(void*));

  nwp->on_edge_change[pos] = callback;
  nwp->on_edge_change_payload[pos] = payload;
  
  nwp->n_on_edge_change++;
}

/*****************
 void DeleteOnWtNetworkEdgeChange

 Delete a specified toggle callback from the list and move the other
 callbacks up the list. Note that both callback and payload pointers
 must match.
*****************/
void DeleteOnWtNetworkEdgeChange(WtNetwork *nwp, OnWtNetworkEdgeChange callback, void *payload){
  unsigned int i;
  for(i = 0; i < nwp->n_on_edge_change; i++)
    if(nwp->on_edge_change[i]==callback && nwp->on_edge_change_payload[i]==payload) break;

  if(i==nwp->n_on_edge_change) error("Attempting to delete a nonexistent callback.");

  memmove(nwp->on_edge_change+i, nwp->on_edge_change+i+1, (nwp->n_on_edge_change-i-1)*sizeof(OnWtNetworkEdgeChange));
  memmove(nwp->on_edge_change_payload+i, nwp->on_edge_change_payload+i+1, (nwp->n_on_edge_change-i-1)*sizeof(void*));

  nwp->n_on_edge_change--;
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
  value of i is EDGECOUNT(nwp).
******************/

/* *** don't forget tail->head, so this function now accepts tail before head */

int WtFindithEdge (Vertex *tail, Vertex *head, double *weight, Edge i, WtNetwork *nwp) {
  Vertex taili=1;
  Edge e;

  /* TODO: This could be speeded up by a factor of 3 or more by starting
     the search from the tail n rather than tail 1 if i > ndyads/2. */

  if (i > EDGECOUNT(nwp) || i<=0)
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
  if(EDGECOUNT(nwp)==0) return(0);
  // FIXME: The constant maxEattempts needs to be tuned.
  const unsigned int maxEattempts=10;
  unsigned int Eattempts = nwp->last_outedge/EDGECOUNT(nwp);
  Edge rane;
  
  if(Eattempts>maxEattempts){
    // If the outedges is too sparse, revert to the old algorithm.
    rane=1 + unif_rand() * EDGECOUNT(nwp);
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
  value of i is (ndyads - EDGECOUNT(nwp)).
******************/

/* This function is not yet written.  It's not clear whether it'll
   be needed. */      
  /* *** but if it is needed, don't forget,  tail -> head */

int WtFindithNonedge (Vertex *tail, Vertex *head, Dyad i, WtNetwork *nwp) {
  Vertex taili=1;
  Edge e;
  Dyad ndyads = DYADCOUNT(nwp);
  
  // If the index is too high or too low, exit immediately.
  if (i > ndyads - EDGECOUNT(nwp) || i<=0)
    return 0;

  /* TODO: This could be speeded up by a factor of 3 or more by starting
     the search from the tail n rather than tail 1 if i > ndyads/2. */

  Vertex nnt;
  while (i > (nnt = nwp->nnodes - (nwp->bipartite ? nwp->bipartite : (nwp->directed_flag?1:taili))
	      - nwp->outdegree[taili])) {   // nnt is the number of nonties incident on taili. Note that when network is undirected, tail<head.
    i -= nnt;
    taili++;
  }

  // Now, our tail is taili.

  /* TODO: This could be speeded up by a factor of 3 or more by starting
     the search from the tree maximum rather than minimum (left over) i > outdegree[taili]. */

 Vertex lhead = (
		  nwp->bipartite ? 
		  nwp->bipartite :
		  (nwp->directed_flag ?
		   taili==1 : taili)
		  );
   e = WtEdgetreeMinimum(nwp->outedges,taili);
  Vertex rhead = nwp->outedges[e].value;
  // Note that rhead-lhead-1-(lhead<taili && taili<rhead) is the number of nonties between two successive ties.
  // the -(lhead<taili && taili<rhead) is because (taili,taili) is not a valid nontie and must be skipped.
  // Note that if taili is an isolate, rhead will be 0.
  while (rhead && i > rhead-lhead-1-(lhead<taili && taili<rhead)) {
    i -= rhead-lhead-1-(lhead<taili && taili<rhead);
    lhead = rhead;
    e = WtEdgetreeSuccessor(nwp->outedges, e);
    // If rhead was the highest-indexed head, then e is now 0.
    if(e) rhead = nwp->outedges[e].value;
    else break; // Note that we don't actually need rhead in the final step.
  }

  // Now, the head we are looking for is (left over) i after lhead.

  *tail = taili;
  *head = lhead + i + (nwp->directed_flag && lhead<taili && lhead+i>=taili); // Skip over the (taili,taili) dyad, if the network is directed.

  return 1;
}

/*****************
  int WtGetRandNonedge

  Select an non-edge in the WtNetwork *nwp at random and update the values
  of tail and head appropriately. Return 1 if successful, 0 otherwise.
******************/

/* *** don't forget tail->head, so this function now accepts tail before head */

int WtGetRandNonedge(Vertex *tail, Vertex *head, WtNetwork *nwp) {
  Dyad ndyads = DYADCOUNT(nwp);
  if(ndyads-EDGECOUNT(nwp)==0) return(0);

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
  unsigned int Eattempts = ndyads/(ndyads-EDGECOUNT(nwp));
  
  if(Eattempts>maxEattempts){
    // If the network is too dense, use the deterministic-time method:
    Dyad rane=1 + unif_rand() * (ndyads-EDGECOUNT(nwp));
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
  return nextedge;
}


/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */


/****************
 Edge WtShuffleEdges

 Randomly permute edges in an list.
****************/
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

/****************
 Edge WtDetShuffleEdges

 Deterministically scramble edges in an list.
****************/
void WtDetShuffleEdges(Vertex *tails, Vertex *heads, double *weights, Edge nedges){
  /* *** don't forget,  tail -> head */
  for(Edge i = nedges; i > 0; i--) {
    Edge j = i/2;
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

/****************
 Edge WtDetUnShuffleEdges

 Reverses WtDetShuffleEdges().
****************/
void WtDetUnShuffleEdges(Vertex *tails, Vertex *heads, double *weights, Edge nedges){
  /* *** don't forget,  tail -> head */
  for(Edge i = 1; i <= nedges; i++) {
    Edge j = i/2;
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

/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */

/*****************
 int WtSetEdge

 Set an weighted edge value: set it to its new weight. Create if it
does not exist, destroy by setting to 0. 
*****************/
void WtSetEdge (Vertex tail, Vertex head, double weight, WtNetwork *nwp) 
{
  ENSURE_TH_ORDER;

  if(weight==0){
    // If the function is to set the edge value to 0, just delete it.
    WtDeleteEdgeFromTrees(tail,head,nwp);
  }else{
    // Find the out-edge
    Edge oe=WtEdgetreeSearch(tail,head,nwp->outedges);
    double w = nwp->outedges[oe].weight;
    if(oe){
      // If it exists AND already has the target weight, do nothing.
      if(w==weight) return;
      else{
        for(unsigned int i = 0; i < nwp->n_on_edge_change; i++) nwp->on_edge_change[i](tail, head, weight, nwp->on_edge_change_payload[i], nwp, w);
	// Find the corresponding in-edge.
	Edge ie=WtEdgetreeSearch(head,tail,nwp->inedges);
	nwp->inedges[ie].weight=nwp->outedges[oe].weight=weight;
      }
    }else{
      // Otherwise, create a new edge with that weight.
      WtAddEdgeToTrees(tail,head,weight,nwp);
    }
  }
}

WtNetwork *Redgelist2WtNetwork(SEXP elR, Rboolean empty){
  Vertex e = length(VECTOR_ELT(elR, 0));
  Vertex *tails = empty ? NULL : (Vertex*) INTEGER(VECTOR_ELT(elR, 0));
  Vertex *heads = empty ? NULL : (Vertex*) INTEGER(VECTOR_ELT(elR, 1));
  double *weights = empty ? NULL : REAL(VECTOR_ELT(elR, 2));
  Vertex n = asInteger(getAttrib(elR, install("n")));
  Rboolean directed = asLogical(getAttrib(elR, install("directed")));
  Vertex bipartite = asInteger(getAttrib(elR, install("bipartite")));
  WtNetwork *nwp = WtNetworkInitialize(tails, heads, weights, e, n, directed, bipartite, FALSE, 0, NULL);
  nwp->eattrname = CHAR(STRING_ELT(getAttrib(elR, R_NamesSymbol),2));
  return nwp;
}


SEXP WtNetwork2Redgelist(WtNetwork *nwp){
  SEXP outl = PROTECT(mkNamed(VECSXP, (const char *[]){".tail", ".head", nwp->eattrname, ""}));
  SEXP newnetworktails = PROTECT(allocVector(INTSXP, EDGECOUNT(nwp)));
  SEXP newnetworkheads = PROTECT(allocVector(INTSXP, EDGECOUNT(nwp)));
  SEXP newnetworkweights = PROTECT(allocVector(REALSXP, EDGECOUNT(nwp)));
  WtEdgeTree2EdgeList((Vertex*)INTEGER(newnetworktails),
                      (Vertex*)INTEGER(newnetworkheads),
                      (double*)REAL(newnetworkweights),
                    nwp,EDGECOUNT(nwp));
  SET_VECTOR_ELT(outl, 0, newnetworktails);
  SET_VECTOR_ELT(outl, 1, newnetworkheads);
  SET_VECTOR_ELT(outl, 2, newnetworkweights);
  UNPROTECT(3);

  SEXP rownames = PROTECT(allocVector(INTSXP, EDGECOUNT(nwp)));
  int *r = INTEGER(rownames);
  for(unsigned int i=1; i<=EDGECOUNT(nwp); i++, r++) *r=i; 
  setAttrib(outl, install("row.names"), rownames);
  UNPROTECT(1);

  setAttrib(outl, install("n"), PROTECT(ScalarInteger(nwp->nnodes)));
  setAttrib(outl, install("directed"), PROTECT(ScalarLogical(nwp->directed_flag)));
  setAttrib(outl, install("bipartite"), PROTECT(ScalarInteger(nwp->bipartite)));
  setAttrib(outl, install("response"), PROTECT(mkChar(nwp->eattrname)));
  UNPROTECT(4);

  SEXP class = PROTECT(mkRStrVec((const char*[]){"tibble_edgelist", "edgelist", "tbl_df", "tbl", "data.frame", NULL}));
  classgets(outl, class);
  UNPROTECT(1);

  UNPROTECT(1);
  return outl;
}
