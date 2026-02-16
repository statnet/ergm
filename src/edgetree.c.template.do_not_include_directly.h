/*  File src/edgetree.c.template.do_not_include_directly.h in package ergm,
 *  part of the Statnet suite of packages for network analysis,
 *  https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */
#include "ergm_Rutil.h"

#include "edgetree_inline_template.do_not_include_directly.h"

/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */


/*******************
 ETYPE(Network) ETYPE(NetworkInitialize)

 Initialize, construct binary tree version of network with weights.  Note
 that the 0th ETYPE(TreeNode) in the array is unused and should
 have all its values set to zero

Note: passing nedges > 0 and tails == heads == NULL is OK: it creates an empty network with at least nedges of space preallocated.
*******************/
/* *** don't forget, tail -> head */

ETYPE(Network) *ETYPE(NetworkInitialize_noLT)(Vertex *tails, Vertex *heads, IFEWT(EWTTYPE *weights,)
			       Edge nedges, Vertex nnodes, Rboolean directed_flag, Vertex bipartite) {
  ETYPE(Network) *nwp = R_Calloc(1, ETYPE(Network));

  IFEWT(nwp->eattrname = NULL);

  nwp->last_inedge = nwp->last_outedge = (Edge)nnodes;
  /* R_Calloc will zero the allocated memory for us, probably a lot
     faster. */
  nwp->outdegree = (Vertex *) R_Calloc((nnodes+1), Vertex);
  nwp->indegree  = (Vertex *) R_Calloc((nnodes+1), Vertex);
  nwp->maxedges = MAX(nedges,1)+nnodes+2; /* Maybe larger than needed? */
  nwp->inedges = (ETYPE(TreeNode) *) R_Calloc(nwp->maxedges, ETYPE(TreeNode));
  nwp->outedges = (ETYPE(TreeNode) *) R_Calloc(nwp->maxedges, ETYPE(TreeNode));

  /*Configure a Network*/
  nwp->nnodes = nnodes;
  EDGECOUNT(nwp) = 0; /* Edges will be added one by one */
  nwp->directed_flag=directed_flag;
  nwp->bipartite=bipartite;

  if(nedges == 0) return nwp;

  ETYPE(DetShuffleEdges)(tails,heads,IFEWT(weights,)nedges); /* shuffle to avoid worst-case performance */

  for(Edge i = 0; i < nedges; i++) {
    Vertex tail=tails[i], head=heads[i];
    IFEWT(EWTTYPE w=weights[i];
          if(w==0) continue);
    if (!directed_flag && tail > head)
      ETYPE(AddEdgeToTrees)(head,tail,IFEWT(w,)nwp); /* Undir edges always have tail < head */
    else
      ETYPE(AddEdgeToTrees)(tail,head,IFEWT(w,)nwp);
  }

  ETYPE(DetUnShuffleEdges)(tails,heads,IFEWT(weights,)nedges); /* Unshuffle edges */

  return nwp;
}

/*******************
 void NetworkDestroy
*******************/
void ETYPE(NetworkDestroy)(ETYPE(Network) *nwp) {
  R_Free(nwp->on_edge_change);
  R_Free(nwp->on_edge_change_payload);
  R_Free(nwp->indegree);
  R_Free(nwp->outdegree);
  R_Free(nwp->inedges);
  R_Free(nwp->outedges);
  R_Free(nwp);
}

/******************
 Network ETYPE(NetworkCopy)
*****************/
ETYPE(Network) *ETYPE(NetworkCopy)(ETYPE(Network) *src){
  ETYPE(Network) *dest = R_Calloc(1, ETYPE(Network));

  Vertex nnodes = dest->nnodes = src->nnodes;
  dest->last_inedge = src->last_inedge;
  dest->last_outedge = src->last_outedge;

  dest->outdegree = (Vertex *) R_Calloc((nnodes+1), Vertex);
  memcpy(dest->outdegree, src->outdegree, (nnodes+1)*sizeof(Vertex));
  dest->indegree = (Vertex *) R_Calloc((nnodes+1), Vertex);
  memcpy(dest->indegree, src->indegree, (nnodes+1)*sizeof(Vertex));

  Vertex maxedges = dest->maxedges = src->maxedges;

  dest->inedges = (ETYPE(TreeNode) *) R_Calloc(maxedges, ETYPE(TreeNode));
  memcpy(dest->inedges, src->inedges, maxedges*sizeof(ETYPE(TreeNode)));
  dest->outedges = (ETYPE(TreeNode) *) R_Calloc(maxedges, ETYPE(TreeNode));
  memcpy(dest->outedges, src->outedges, maxedges*sizeof(ETYPE(TreeNode)));

  dest->directed_flag = src->directed_flag;
  dest->bipartite = src->bipartite;

  EDGECOUNT(dest) = EDGECOUNT(src);

  return dest;
}

/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */



/*****************
 int ETYPE(ToggleEdge)

 Toggle an edge:  Set it to the opposite of its current
 value.  Return 1 if edge added, 0 if deleted.
*****************/

/* *** don't forget tail->head, so this function now accepts tail before head */

int ETYPE(ToggleEdge) (Vertex tail, Vertex head, IFEWT(EWTTYPE weight,) ETYPE(Network) *nwp)
{
  /* don't forget tails < heads now for undirected networks */
  ENSURE_TH_ORDER;
  if(ETYPE(DeleteEdgeFromTrees)(tail,head,nwp))
    return 0;
  else{
    ETYPE(AddEdgeToTrees)(tail,head,IFEWT(weight,)nwp);
    return 1;
  }
}

/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */

/*****************
 void AddEdgeToTrees

 Add an edge from tail to head; note that the function assumes that it
 is legal and therefore does not have a return value. Since each
 "edge" should be added to both the list of outedges and the list of
 inedges, this actually involves two calls to AddHalfedgeToTree (hence
 "Trees" instead of "Tree" in the name of this function).
*****************/

/* *** don't forget tail->head, so this function now accepts tail before head */

void ETYPE(AddEdgeToTrees)(Vertex tail, Vertex head, IFEWT(EWTTYPE weight,) ETYPE(Network) *nwp){
#ifndef NDEBUG
  if(!nwp->directed_flag && tail>head)
    error("ETYPE(AddEdgeToTrees)() called for an undirected network with tail>head. Note that this produces an error only if compiling with NDEBUG macro unset and silently produces undefined behavior otherwise.");
  if(ETYPE(EdgetreeSearch)(tail, head, nwp->outedges)||ETYPE(EdgetreeSearch)(head, tail, nwp->inedges)) error("ETYPE(AddEdgeToTrees)() called for an extant edge. Note that this produces an error only if compiling with NDEBUG macro unset and silently produces undefined behavior otherwise.");
#endif // NDEBUG
  for(unsigned int i = 0; i < nwp->n_on_edge_change; i++) nwp->on_edge_change[i](tail, head, IFEWT(weight,) nwp->on_edge_change_payload[i], nwp, 0);

  ETYPE(AddHalfedgeToTree)(tail, head, IFEWT(weight,) nwp->outedges, &(nwp->last_outedge));
  ETYPE(AddHalfedgeToTree)(head, tail, IFEWT(weight,) nwp->inedges, &(nwp->last_inedge));
  ++nwp->outdegree[tail];
  ++nwp->indegree[head];
  ++EDGECOUNT(nwp);
  ETYPE(CheckEdgetreeFull)(nwp);
}

/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */

/*****************
 int ETYPE(DeleteEdgeFromTrees)

 Find and delete the edge from tail to head.
 Return 1 if successful, 0 otherwise.  As with AddEdgeToTrees, this must
 be done once for outedges and once for inedges.
*****************/

/* *** don't forget tail->head, so this function now accepts tail before head */

int ETYPE(DeleteEdgeFromTrees)(Vertex tail, Vertex head, ETYPE(Network) *nwp){
  Edge zth, zht;
  if((zth=ETYPE(EdgetreeSearch)(tail, head, nwp->outedges))&&(zht=ETYPE(EdgetreeSearch)(head, tail, nwp->inedges))){
    if(nwp->n_on_edge_change){
      IFEWT(EWTTYPE w = nwp->outedges[zth].weight);
      for(unsigned int i = 0; i < nwp->n_on_edge_change; i++) nwp->on_edge_change[i](tail, head, IFEWT(0,) nwp->on_edge_change_payload[i], nwp, IFELSEEWT(w,TRUE));
    }
    ETYPE(DeleteHalfedgeFromTreeAt)(tail, head, nwp->outedges,&(nwp->last_outedge), zth);
    ETYPE(DeleteHalfedgeFromTreeAt)(head, tail, nwp->inedges, &(nwp->last_inedge), zht);
    --nwp->outdegree[tail];
    --nwp->indegree[head];
    --EDGECOUNT(nwp);
    return 1;
  }
  return 0;
}

/*****************
 void ETYPE(AddOn, NetworkEdgeChange)

 Insert a specified toggle callback at the specified position.
*****************/
void ETYPE(AddOn, NetworkEdgeChange)(ETYPE(Network) *nwp, ETYPE(On, NetworkEdgeChange) callback, void *payload, unsigned int pos){
  if(nwp->n_on_edge_change+1 > nwp->max_on_edge_change){
    nwp->max_on_edge_change = MAX(nwp->max_on_edge_change,1)*2;
    nwp->on_edge_change = R_Realloc(nwp->on_edge_change, nwp->max_on_edge_change, ETYPE(On, NetworkEdgeChange));
    nwp->on_edge_change_payload = R_Realloc(nwp->on_edge_change_payload, nwp->max_on_edge_change, void*);
  }

  pos = MIN(pos, nwp->n_on_edge_change); // Last position.
  // Move everything down the list.
  memmove(nwp->on_edge_change+pos+1, nwp->on_edge_change+pos, (nwp->n_on_edge_change-pos)*sizeof(ETYPE(On, NetworkEdgeChange)));
  memmove(nwp->on_edge_change_payload+pos+1, nwp->on_edge_change_payload+pos, (nwp->n_on_edge_change-pos)*sizeof(void*));

  nwp->on_edge_change[pos] = callback;
  nwp->on_edge_change_payload[pos] = payload;

  nwp->n_on_edge_change++;
}

/*****************
 void ETYPE(DeleteOn, NetworkEdgeChange)

 Delete a specified toggle callback from the list and move the other
 callbacks up the list. Note that both callback and payload pointers
 must match.
*****************/
void ETYPE(DeleteOn, NetworkEdgeChange)(ETYPE(Network) *nwp, ETYPE(On, NetworkEdgeChange) callback, void *payload){
  unsigned int i;
  for(i = 0; i < nwp->n_on_edge_change; i++)
    if(nwp->on_edge_change[i]==callback && nwp->on_edge_change_payload[i]==payload) break;

  if(i==nwp->n_on_edge_change) error("Attempting to delete a nonexistent callback.");

  memmove(nwp->on_edge_change+i, nwp->on_edge_change+i+1, (nwp->n_on_edge_change-i-1)*sizeof(ETYPE(On, NetworkEdgeChange)));
  memmove(nwp->on_edge_change_payload+i, nwp->on_edge_change_payload+i+1, (nwp->n_on_edge_change-i-1)*sizeof(void*));

  nwp->n_on_edge_change--;
}

/*****************
 void ETYPE(NetworkEdgeList)

 Print an edgelist for a network
*****************/
void ETYPE(NetworkEdgeList)(ETYPE(Network) *nwp) {
  Vertex i;
  for(i=1; i<=nwp->nnodes; i++) {
    Rprintf("Node %d:\n  ", i);
    ETYPE(InOrderTreeWalk)(nwp->outedges, i);
    Rprintf("\n");
  }
}

/*****************
 void ETYPE(InOrderTreeWalk)

 Diagnostic routine that prints the nodes in the tree rooted
 at edges[x], in increasing order of their values.
*****************/
void ETYPE(InOrderTreeWalk)(ETYPE(TreeNode) *edges, Edge x) {
  if (x != 0) {
    ETYPE(InOrderTreeWalk)(edges, (edges+x)->left);
    IFELSEEWT(Rprintf(" %d:%f ",(edges+x)->value, (edges+x)->weight),
              Rprintf(" %d ",(edges+x)->value));
    ETYPE(InOrderTreeWalk)(edges, (edges+x)->right);
  }
}

/*****************
  int ETYPE(FindithEdge)

  Find the ith edge in the ETYPE(Network) *nwp and update the values of
  tail, head, and weight appropriately. If the value passed to tail,
  head, or weight is NULL, it is not updated, so it is possible to
  only obtain what is needed. Return 1 if successful, 0 otherwise.
  Note that i is numbered from 1, not 0.  Thus, the maximum possible
  value of i is EDGECOUNT(nwp).
******************/

/* *** don't forget tail->head, so this function now accepts tail before head */

int ETYPE(FindithEdge) (Vertex *tail, Vertex *head, IFEWT(EWTTYPE *weight,) Edge i, ETYPE(Network) *nwp) {
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

  e=ETYPE(EdgetreeMinimum)(nwp->outedges,taili);
  while (i-- > 1) {
    e=ETYPE(EdgetreeSuccessor)(nwp->outedges, e);
  }
  if(tail) *tail = taili;
  if(head) *head = nwp->outedges[e].value;
  IFEWT(if(weight) *weight = nwp->outedges[e].weight);
  return 1;
}

/*****************
  int GetRandEdge

  Select an edge in the Network *nwp at random and update the values
  of tail and head appropriately. Return 1 if successful, 0 otherwise.
******************/

/* *** don't forget tail->head, so this function now accepts tail before head */

int ETYPE(GetRandEdge)(Vertex *tail, Vertex *head, IFEWT(EWTTYPE *weight,) ETYPE(Network) *nwp) {
  if(EDGECOUNT(nwp)==0) return(0);
  // FIXME: The constant maxEattempts needs to be tuned.
  const unsigned int maxEattempts=10;
  unsigned int Eattempts = nwp->last_outedge/EDGECOUNT(nwp);
  Edge rane;

  if(Eattempts>maxEattempts){
    // If the outedges is too sparse, revert to the old algorithm.
    rane=1 + unif_rand() * EDGECOUNT(nwp);
    ETYPE(FindithEdge)(tail, head, IFEWT(weight,) rane, nwp);
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
    IFEWT(if(weight) *weight=nwp->outedges[rane].weight);

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
  int ETYPE(FindithNonedge)

  Find the ith nonedge in the ETYPE(Network) *nwp and
  update the values of tail and head appropriately.  Return
  1 if successful, 0 otherwise.
  Note that i is numbered from 1, not 0.  Thus, the maximum possible
  value of i is (ndyads - EDGECOUNT(nwp)).
******************/

/* This function is not yet written.  It's not clear whether it'll
   be needed. */
  /* *** but if it is needed, don't forget,  tail -> head */

int ETYPE(FindithNonedge) (Vertex *tail, Vertex *head, Dyad i, ETYPE(Network) *nwp) {
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

  // If taili 1, then head cannot be 1. If undirected, the smallest it can be is taili+1. If bipartite, the smallest it can be is nwp->bipartite+1.
  Vertex lhead = (
		  nwp->bipartite ?
		  nwp->bipartite :
		  (nwp->directed_flag ?
		   taili==1 : taili)
		  );
   e = ETYPE(EdgetreeMinimum)(nwp->outedges,taili);
  Vertex rhead = nwp->outedges[e].value;
  // Note that rhead-lhead-1-(lhead<taili && taili<rhead) is the number of nonties between two successive ties.
  // the -(lhead<taili && taili<rhead) is because (taili,taili) is not a valid nontie and must be skipped.
  // Note that if taili is an isolate, rhead will be 0.
  while (rhead && i > rhead-lhead-1-(lhead<taili && taili<rhead)) {
    i -= rhead-lhead-1-(lhead<taili && taili<rhead);
    lhead = rhead;
    e = ETYPE(EdgetreeSuccessor)(nwp->outedges, e);
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
  int ETYPE(GetRandNonedge)

  Select an non-edge in the ETYPE(Network) *nwp at random and update the values
  of tail and head appropriately. Return 1 if successful, 0 otherwise.
******************/

/* *** don't forget tail->head, so this function now accepts tail before head */

int ETYPE(GetRandNonedge)(Vertex *tail, Vertex *head, ETYPE(Network) *nwp) {
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
    ETYPE(FindithNonedge)(tail, head, rane, nwp);
  }else{
    do{
      GetRandDyad(tail, head, nwp);
    }while(ETYPE(EdgetreeSearch)(*tail, *head, nwp->outedges));
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
Edge ETYPE(EdgeTree2EdgeList)(Vertex *tails, Vertex *heads, IFEWT(EWTTYPE *weights,) ETYPE(Network) *nwp, Edge nmax){
  Edge nextedge=0;

  /* *** don't forget,  tail -> head */
  for (Vertex v=1; v<=nwp->nnodes; v++){
    for(Vertex e = ETYPE(EdgetreeMinimum)(nwp->outedges, v);
	nwp->outedges[e].value != 0 && nextedge < nmax;
	e = ETYPE(EdgetreeSuccessor)(nwp->outedges, e)){
      tails[nextedge] = v;
      heads[nextedge] = nwp->outedges[e].value;
      IFEWT(if(weights) weights[nextedge] = nwp->outedges[e].weight);
      nextedge++;
    }
  }
  return nextedge;
}


/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */


/****************
 Edge ETYPE(ShuffleEdges)

 Randomly permute edges in an list.
****************/
void ETYPE(ShuffleEdges)(Vertex *tails, Vertex *heads, IFEWT(EWTTYPE *weights,) Edge nedges){
  /* *** don't forget,  tail -> head */
  for(Edge i = nedges; i > 0; i--) {
    Edge j = i * unif_rand();
    Vertex tail = tails[j];
    Vertex head = heads[j];
    IFEWT(EWTTYPE w = weights[j]);
    tails[j] = tails[i-1];
    heads[j] = heads[i-1];
    IFEWT(weights[j] = weights[i-1]);
    tails[i-1] = tail;
    heads[i-1] = head;
    IFEWT(weights[i-1] = w);
  }
}

/****************
 Edge ETYPE(DetShuffleEdges)

 Deterministically scramble edges in an list.
****************/
void ETYPE(DetShuffleEdges)(Vertex *tails, Vertex *heads, IFEWT(EWTTYPE *weights,) Edge nedges){
  /* *** don't forget,  tail -> head */
  for(Edge i = nedges; i > 0; i--) {
    Edge j = i/2;
    Vertex tail = tails[j];
    Vertex head = heads[j];
    IFEWT(EWTTYPE w = weights[j]);
    tails[j] = tails[i-1];
    heads[j] = heads[i-1];
    IFEWT(weights[j] = weights[i-1]);
    tails[i-1] = tail;
    heads[i-1] = head;
    IFEWT(weights[i-1] = w);
  }
}

/****************
 Edge ETYPE(DetUnShuffleEdges)

 Reverses ETYPE(DetShuffleEdges)().
****************/
void ETYPE(DetUnShuffleEdges)(Vertex *tails, Vertex *heads, IFEWT(EWTTYPE *weights,) Edge nedges){
  /* *** don't forget,  tail -> head */
  for(Edge i = 1; i <= nedges; i++) {
    Edge j = i/2;
    Vertex tail = tails[j];
    Vertex head = heads[j];
    IFEWT(EWTTYPE w = weights[j]);
    tails[j] = tails[i-1];
    heads[j] = heads[i-1];
    IFEWT(weights[j] = weights[i-1]);
    tails[i-1] = tail;
    heads[i-1] = head;
    IFEWT(weights[i-1] = w);
  }
}

/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */

/*****************
 int ETYPE(SetEdge)

 Set an weighted edge value: set it to its new weight. Create if it
does not exist, destroy by setting to 0.
*****************/
void ETYPE(SetEdge) (Vertex tail, Vertex head, EWTTYPE weight, ETYPE(Network) *nwp)
{
  ENSURE_TH_ORDER;

  if(weight==0){
    // If the function is to set the edge value to 0, just delete it.
    ETYPE(DeleteEdgeFromTrees)(tail,head,nwp);
  }else{
    // Find the out-edge
    Edge oe=ETYPE(EdgetreeSearch)(tail,head,nwp->outedges);
    IFEWT(EWTTYPE w = nwp->outedges[oe].weight);
    if(oe == 0) ETYPE(AddEdgeToTrees)(tail,head,IFEWT(weight,)nwp);
    IFEWT(else{
      // If it exists AND already has the target weight, do nothing.
      if(w==weight) return;
      else{
        for(unsigned int i = 0; i < nwp->n_on_edge_change; i++) nwp->on_edge_change[i](tail, head, weight, nwp->on_edge_change_payload[i], nwp, w);
	// Find the corresponding in-edge.
	Edge ie=ETYPE(EdgetreeSearch)(head,tail,nwp->inedges);
	nwp->inedges[ie].weight=nwp->outedges[oe].weight=weight;
      }
    })
  }
}

ETYPE(Network) *ETYPE(Redgelist2, Network)(SEXP elR, Rboolean empty){
  Vertex e = empty ? 0 : length(VECTOR_ELT(elR, 0));
  Vertex *tails = empty ? NULL : (Vertex*) INTEGER(VECTOR_ELT(elR, 0));
  Vertex *heads = empty ? NULL : (Vertex*) INTEGER(VECTOR_ELT(elR, 1));
  IFEWT(EWTTYPE *weights = empty ? NULL : REAL(VECTOR_ELT(elR, 2)));
  Vertex n = asInteger(getAttrib(elR, install("n")));
  Rboolean directed = asLogical(getAttrib(elR, install("directed")));
  Vertex bipartite = asInteger(getAttrib(elR, install("bipartite")));
  ETYPE(Network) *nwp = ETYPE(NetworkInitialize)(tails, heads, IFEWT(weights,) e, n, directed, bipartite);
  IFEWT(nwp->eattrname = CHAR(STRING_ELT(getAttrib(elR, R_NamesSymbol),2)));
  return nwp;
}


SEXP ETYPE(Network2Redgelist)(ETYPE(Network) *nwp){
  SEXP outl = PROTECT(mkNamed(VECSXP, (const char *[]){".tail", ".head", IFEWT(nwp->eattrname,) ""}));
  SEXP newnetworktails = PROTECT(allocVector(INTSXP, EDGECOUNT(nwp)));
  SEXP newnetworkheads = PROTECT(allocVector(INTSXP, EDGECOUNT(nwp)));
  IFEWT(SEXP newnetworkweights = PROTECT(allocVector(REALSXP, EDGECOUNT(nwp))));
  ETYPE(EdgeTree2EdgeList)((Vertex*)INTEGER(newnetworktails),
                      (Vertex*)INTEGER(newnetworkheads),
                      IFEWT((EWTTYPE*)REAL(newnetworkweights),)
                    nwp,EDGECOUNT(nwp));
  SET_VECTOR_ELT(outl, 0, newnetworktails);
  SET_VECTOR_ELT(outl, 1, newnetworkheads);
  IFEWT(SET_VECTOR_ELT(outl, 2, newnetworkweights));
  UNPROTECT(IFELSEEWT(3,2));

  SEXP rownames = PROTECT(allocVector(INTSXP, EDGECOUNT(nwp)));
  int *r = INTEGER(rownames);
  for(unsigned int i=1; i<=EDGECOUNT(nwp); i++, r++) *r=i;
  setAttrib(outl, install("row.names"), rownames);
  UNPROTECT(1);

  setAttrib(outl, install("n"), PROTECT(ScalarInteger(nwp->nnodes)));
  setAttrib(outl, install("directed"), PROTECT(ScalarLogical(nwp->directed_flag)));
  setAttrib(outl, install("bipartite"), PROTECT(ScalarInteger(nwp->bipartite)));
  IFEWT(setAttrib(outl, install("response"), PROTECT(mkChar(nwp->eattrname))));
  UNPROTECT(IFELSEEWT(4,3));

  SEXP class = PROTECT(mkRStrVec((const char*[]){"tibble_edgelist", "edgelist", "tbl_df", "tbl", "data.frame", NULL}));
  classgets(outl, class);
  UNPROTECT(1);

  UNPROTECT(1);
  return outl;
}
