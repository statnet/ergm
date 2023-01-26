/*  File src/edgetree.c in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2023 Statnet Commons
 */
#include "ergm_edgetree.h"
#include "ergm_Rutil.h"

#include "edgetree_inline.do_not_include_directly.h"

/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */


/*******************
 Network NetworkInitialize

 Initialize, construct binary tree version of network.  Note
 that the 0th TreeNode in the array is unused and should 
 have all its values set to zero

Note: passing nedges > 0 and tails == heads == weights == NULL is OK: it creates an empty network with at least nedges of space preallocated.
*******************/
/* *** don't forget, tail -> head */

Network *NetworkInitialize(Vertex *tails, Vertex *heads, Edge nedges, 
			   Vertex nnodes, int directed_flag, Vertex bipartite,
			   int lasttoggle_flag, int time, int *lasttoggle) {
  
  Network *nwp = Calloc(1, Network);

  nwp->last_inedge = nwp->last_outedge = (Edge)nnodes;
  /* Calloc will zero the allocated memory for us, probably a lot
     faster. */
  nwp->outdegree = (Vertex *) Calloc((nnodes+1), Vertex);
  nwp->indegree  = (Vertex *) Calloc((nnodes+1), Vertex);
  nwp->maxedges = MAX(nedges,1)+nnodes+2; /* Maybe larger than needed? */
  nwp->inedges = (TreeNode *) Calloc(nwp->maxedges, TreeNode);
  nwp->outedges = (TreeNode *) Calloc(nwp->maxedges, TreeNode);

  if(lasttoggle_flag) error("The lasttoggle API has been removed from ergm.");

  /*Configure a Network*/
  nwp->nnodes = nnodes;
  EDGECOUNT(nwp) = 0; /* Edges will be added one by one */
  nwp->directed_flag=directed_flag;
  nwp->bipartite=bipartite;

  if(tails==NULL && heads==NULL) return nwp;

  DetShuffleEdges(tails,heads,nedges); /* shuffle to avoid worst-case performance */

  for(Edge i = 0; i < nedges; i++) {
    Vertex tail=tails[i], head=heads[i];
    if (!directed_flag && tail > head) 
      AddEdgeToTrees(head,tail,nwp); /* Undir edges always have tail < head */ 
    else 
      AddEdgeToTrees(tail,head,nwp);
  }

  DetUnShuffleEdges(tails,heads,nedges); /* Unshuffle edges */

  return nwp;
}


/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */

/*Takes vectors of doubles for edges; used only when constructing from inputparams. */
Network *NetworkInitializeD(double *tails, double *heads, Edge nedges,
			  Vertex nnodes, int directed_flag, Vertex bipartite,
			  int lasttoggle_flag, int time, int *lasttoggle) {

  /* *** don't forget, tail -> head */

  Vertex *itails=(Vertex*)Calloc(nedges, Vertex);
  Vertex *iheads=(Vertex*)Calloc(nedges, Vertex);
  
  for(Edge i=0; i<nedges; i++){
    itails[i]=tails[i];
    iheads[i]=heads[i];
  }

  Network *nwp=NetworkInitialize(itails,iheads,nedges,nnodes,directed_flag,bipartite,lasttoggle_flag, time, lasttoggle);

  Free(itails);
  Free(iheads);
  return nwp;
}

/*******************
 void NetworkDestroy
*******************/
void NetworkDestroy(Network *nwp) {
  Free(nwp->on_edge_change);
  Free(nwp->on_edge_change_payload);
  Free(nwp->indegree);
  Free(nwp->outdegree);
  Free(nwp->inedges);
  Free(nwp->outedges);
  Free(nwp);
}

/******************
 Network NetworkCopy
*****************/
Network *NetworkCopy(Network *src){
  Network *dest = Calloc(1, Network);
  
  Vertex nnodes = dest->nnodes = src->nnodes;
  dest->last_inedge = src->last_inedge;
  dest->last_outedge = src->last_outedge;

  dest->outdegree = (Vertex *) Calloc((nnodes+1), Vertex);
  memcpy(dest->outdegree, src->outdegree, (nnodes+1)*sizeof(Vertex));
  dest->indegree = (Vertex *) Calloc((nnodes+1), Vertex);
  memcpy(dest->indegree, src->indegree, (nnodes+1)*sizeof(Vertex));

  Vertex maxedges = dest->maxedges = src->maxedges;

  dest->inedges = (TreeNode *) Calloc(maxedges, TreeNode);
  memcpy(dest->inedges, src->inedges, maxedges*sizeof(TreeNode));
  dest->outedges = (TreeNode *) Calloc(maxedges, TreeNode);
  memcpy(dest->outedges, src->outedges, maxedges*sizeof(TreeNode));

  dest->directed_flag = src->directed_flag;
  dest->bipartite = src->bipartite;

  EDGECOUNT(dest) = EDGECOUNT(src);

  return dest;
}

/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */



/*****************
 Edge ToggleEdge

 Toggle an edge:  Set it to the opposite of its current
 value.  Return 1 if edge added, 0 if deleted.
*****************/

/* *** don't forget tail->head, so this function now accepts tail before head */

int ToggleEdge (Vertex tail, Vertex head, Network *nwp) 
{
  /* don't forget tails < heads now for undirected networks */
  ENSURE_TH_ORDER;
  if(DeleteEdgeFromTrees(tail,head,nwp))
    return 0;
  else{
    AddEdgeToTrees(tail,head,nwp);
    return 1;
  }
}


/*****************
 Edge ToggleKnownEdge

 Toggle an edge whose status is known:  Set it to the opposite of its current
 value.  Return 1 if edge added, 0 if deleted.
*****************/

/* *** don't forget tail->head, so this function now accepts tail before head */

void ToggleKnownEdge (Vertex tail, Vertex head, Network *nwp, Rboolean edgestate)
{
#ifdef DEBUG
  if((EdgetreeSearch(tail, head, nwp->outedges)!=0) != edgestate)
    error("ToggleKnownEdge() called with an incorrect edgestate. Note that this produces an error only if compiling with DEBUG macro set and silently produces undefined behavior otherwise.");
#endif // DEBUG
  ENSURE_TH_ORDER;
  if (edgestate){
    DeleteEdgeFromTrees(tail,head,nwp);
  }else{
    AddEdgeToTrees(tail,head,nwp);
  }
}


/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */


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

void AddEdgeToTrees(Vertex tail, Vertex head, Network *nwp){
#ifdef DEBUG
  if(!nwp->directed_flag && tail>head)
    error("AddEdgeToTrees() called for an undirected network with tail>head. Note that this produces an error only if compiling with DEBUG macro set and silently produces undefined behavior otherwise.");
  if(EdgetreeSearch(tail, head, nwp->outedges)||EdgetreeSearch(head, tail, nwp->inedges))
    error("AddEdgeToTrees() called for an extant edge. Note that this produces an error only if compiling with DEBUG macro set and silently produces undefined behavior otherwise.");
#endif // DEBUG
  for(unsigned int i = 0; i < nwp->n_on_edge_change; i++) nwp->on_edge_change[i](tail, head, nwp->on_edge_change_payload[i], nwp, FALSE);

  AddHalfedgeToTree(tail, head, nwp->outedges, &(nwp->last_outedge));
  AddHalfedgeToTree(head, tail, nwp->inedges, &(nwp->last_inedge));
  ++nwp->outdegree[tail];
  ++nwp->indegree[head];
  ++EDGECOUNT(nwp);
  CheckEdgetreeFull(nwp);
}

/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */

/*****************
 int DeleteEdgeFromTrees

 Find and delete the edge from tail to head.  
 Return 1 if successful, 0 otherwise.  As with AddEdgeToTrees, this must
 be done once for outedges and once for inedges.
*****************/

/* *** don't forget tail->head, so this function now accepts tail before head */

int DeleteEdgeFromTrees(Vertex tail, Vertex head, Network *nwp){
  Edge zth, zht;
  if((zth=EdgetreeSearch(tail, head, nwp->outedges))&&(zht=EdgetreeSearch(head, tail, nwp->inedges))){
    if(nwp->n_on_edge_change){
      for(unsigned int i = 0; i < nwp->n_on_edge_change; i++) nwp->on_edge_change[i](tail, head, nwp->on_edge_change_payload[i], nwp, TRUE);
    }
    DeleteHalfedgeFromTreeAt(tail, head, nwp->outedges,&(nwp->last_outedge), zth);
    DeleteHalfedgeFromTreeAt(head, tail, nwp->inedges, &(nwp->last_inedge), zht);
    --nwp->outdegree[tail];
    --nwp->indegree[head];
    --EDGECOUNT(nwp);
    return 1;
  }
  return 0;
}


/*****************
 void AddOnNetworkEdgeChange

 Insert a specified toggle callback at the specified position.
*****************/
void AddOnNetworkEdgeChange(Network *nwp, OnNetworkEdgeChange callback, void *payload, unsigned int pos){
  if(nwp->n_on_edge_change+1 > nwp->max_on_edge_change){
    nwp->max_on_edge_change = MAX(nwp->max_on_edge_change,1)*2;
    nwp->on_edge_change = Realloc(nwp->on_edge_change, nwp->max_on_edge_change, OnNetworkEdgeChange);
    nwp->on_edge_change_payload = Realloc(nwp->on_edge_change_payload, nwp->max_on_edge_change, void*);
  }

  pos = MIN(pos, nwp->n_on_edge_change); // Last position.
  // Move everything down the list.
  memmove(nwp->on_edge_change+pos+1, nwp->on_edge_change+pos, (nwp->n_on_edge_change-pos)*sizeof(OnNetworkEdgeChange));
  memmove(nwp->on_edge_change_payload+pos+1, nwp->on_edge_change_payload+pos, (nwp->n_on_edge_change-pos)*sizeof(void*));

  nwp->on_edge_change[pos] = callback;
  nwp->on_edge_change_payload[pos] = payload;
  
  nwp->n_on_edge_change++;
}

/*****************
 void DeleteOnNetworkEdgeChange

 Delete a specified toggle callback from the list and move the other
 callbacks up the list. Note that both callback and payload pointers
 must match.
*****************/
void DeleteOnNetworkEdgeChange(Network *nwp, OnNetworkEdgeChange callback, void *payload){
  unsigned int i;
  for(i = 0; i < nwp->n_on_edge_change; i++)
    if(nwp->on_edge_change[i]==callback && nwp->on_edge_change_payload[i]==payload) break;

  if(i==nwp->n_on_edge_change) error("Attempting to delete a nonexistent callback.");

  memmove(nwp->on_edge_change+i, nwp->on_edge_change+i+1, (nwp->n_on_edge_change-i-1)*sizeof(OnNetworkEdgeChange));
  memmove(nwp->on_edge_change_payload+i, nwp->on_edge_change_payload+i+1, (nwp->n_on_edge_change-i-1)*sizeof(void*));

  nwp->n_on_edge_change--;
}


/*****************
 void printedge

 Diagnostic routine that prints out the contents
 of the specified TreeNode (used for debugging).  
*****************/
void printedge(Edge e, TreeNode *edges){
  Rprintf("Edge structure [%d]:\n",e);
  Rprintf("\t.value=%d\n",edges[e].value);
  Rprintf("\t.parent=%d\n",edges[e].parent);
  Rprintf("\t.left=%d\n",edges[e].left);
  Rprintf("\t.right=%d\n",edges[e].right);
}

/*****************
 void NetworkEdgeList

 Print an edgelist for a network
*****************/
void NetworkEdgeList(Network *nwp) {
  Vertex i;
  for(i=1; i<=nwp->nnodes; i++) {
    Rprintf("Node %d:\n  ", i);
    InOrderTreeWalk(nwp->outedges, i);
    Rprintf("\n");
  }
}

/*****************
 void InOrderTreeWalk

 Diagnostic routine that prints the nodes in the tree rooted
 at edges[x], in increasing order of their values.
*****************/
void InOrderTreeWalk(TreeNode *edges, Edge x) {
  if (x != 0) {
    InOrderTreeWalk(edges, (edges+x)->left);
    /*    printedge(x, edges); */
    Rprintf(" %d ",(edges+x)->value); 
    InOrderTreeWalk(edges, (edges+x)->right);
  }
}

/*****************
 Edge DesignMissing (see EdgetreeSearch)
*****************/
Edge DesignMissing (Vertex a, Vertex b, Network *mnwp) {
  int miss;
  miss = EdgetreeSearch(a,b,mnwp->outedges);
  if(mnwp->directed_flag){
    miss += EdgetreeSearch(a,b,mnwp->inedges);
  }
  return(miss);
}


/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */


/*****************
  int FindithEdge

  Find the ith edge in the Network *nwp and
  update the values of tail and head appropriately.  Return
  1 if successful, 0 otherwise.  
  Note that i is numbered from 1, not 0.  Thus, the maximum possible
  value of i is EDGECOUNT(nwp).
******************/

/* *** don't forget tail->head, so this function now accepts tail before head */

int FindithEdge (Vertex *tail, Vertex *head, Edge i, Network *nwp) {
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

  e=EdgetreeMinimum(nwp->outedges,taili);
  while (i-- > 1) {
    e=EdgetreeSuccessor(nwp->outedges, e);
  }
  *tail = taili;
  *head = nwp->outedges[e].value;
  return 1;
}

/*****************
  int GetRandEdge

  Select an edge in the Network *nwp at random and update the values
  of tail and head appropriately. Return 1 if successful, 0 otherwise.
******************/

/* *** don't forget tail->head, so this function now accepts tail before head */

int GetRandEdge(Vertex *tail, Vertex *head, Network *nwp) {
  if(EDGECOUNT(nwp)==0) return(0);
  // FIXME: The constant maxEattempts needs to be tuned.
  const unsigned int maxEattempts=10;
  unsigned int Eattempts = nwp->last_outedge/EDGECOUNT(nwp);
  Edge rane;
  
  if(Eattempts>maxEattempts){
    // If the outedges is too sparse, revert to the old algorithm.
    rane=1 + unif_rand() * EDGECOUNT(nwp);
    FindithEdge(tail, head, rane, nwp);
  }else{
    // Otherwise, find a TreeNode which has a head.
    do{
      // Note that the outedges array has maxedges elements, but the
      // 0th one is always blank, and those with index >
      // nwp->last_outedge are blank as well, so we need to generate
      // an index from 1 through nwp->last_outedge (inclusive).
      rane = 1 + unif_rand() * nwp->last_outedge;
    }while((*head=nwp->outedges[rane].value)==0); // Form the head, while we are at it.
    
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
  int FindithNonedge

  Find the ith nonedge in the Network *nwp and
  update the values of tail and head appropriately.  Return
  1 if successful, 0 otherwise.  
  Note that i is numbered from 1, not 0.  Thus, the maximum possible
  value of i is (ndyads - EDGECOUNT(nwp)).
******************/

  /* *** don't forget,  tail -> head */

int FindithNonedge (Vertex *tail, Vertex *head, Dyad i, Network *nwp) {
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
  e = EdgetreeMinimum(nwp->outedges,taili);
  Vertex rhead = nwp->outedges[e].value;
  // Note that rhead-lhead-1-(lhead<taili && taili<rhead) is the number of nonties between two successive ties.
  // the -(lhead<taili && taili<rhead) is because (taili,taili) is not a valid nontie and must be skipped.
  // Note that if taili is an isolate, rhead will be 0.
  while (rhead && i > rhead-lhead-1-(lhead<taili && taili<rhead)) {
    i -= rhead-lhead-1-(lhead<taili && taili<rhead);
    lhead = rhead;
    e = EdgetreeSuccessor(nwp->outedges, e);
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
  int GetRandNonedge

  Select an non-edge in the Network *nwp at random and update the values
  of tail and head appropriately. Return 1 if successful, 0 otherwise.
******************/

/* *** don't forget tail->head, so this function now accepts tail before head */

int GetRandNonedge(Vertex *tail, Vertex *head, Network *nwp) {
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
    FindithNonedge(tail, head, rane, nwp);
  }else{
    do{
      GetRandDyad(tail, head, nwp);
    }while(EdgetreeSearch(*tail, *head, nwp->outedges));
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
Edge EdgeTree2EdgeList(Vertex *tails, Vertex *heads, Network *nwp, Edge nmax){
  Edge nextedge=0;

  /* *** don't forget,  tail -> head */
  for (Vertex v=1; v<=nwp->nnodes; v++){
    for(Vertex e = EdgetreeMinimum(nwp->outedges, v);
	nwp->outedges[e].value != 0 && nextedge < nmax;
	e = EdgetreeSuccessor(nwp->outedges, e)){
      tails[nextedge] = v;
      heads[nextedge] = nwp->outedges[e].value;
      nextedge++;
    }
  }
  return nextedge;
}


/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */


/****************
 Edge ShuffleEdges

 Randomly permute edges in an list.
****************/
void ShuffleEdges(Vertex *tails, Vertex *heads, Edge nedges){
  /* *** don't forget,  tail -> head */
  for(Edge i = nedges; i > 0; i--) {
    Edge j = i * unif_rand();
    Vertex tail = tails[j];
    Vertex head = heads[j];
    tails[j] = tails[i-1];
    heads[j] = heads[i-1];
    tails[i-1] = tail;
    heads[i-1] = head;
  }
}

/****************
 Edge DetShuffleEdges

 Deterministically scramble edges in an list.
****************/
void DetShuffleEdges(Vertex *tails, Vertex *heads, Edge nedges){
  /* *** don't forget,  tail -> head */
  for(Edge i = nedges; i > 0; i--) {
    Edge j = i/2;
    Vertex tail = tails[j];
    Vertex head = heads[j];
    tails[j] = tails[i-1];
    heads[j] = heads[i-1];
    tails[i-1] = tail;
    heads[i-1] = head;
  }
}

/****************
 Edge DetUnShuffleEdges

 Reverses DetShuffleEdges().
****************/
void DetUnShuffleEdges(Vertex *tails, Vertex *heads, Edge nedges){
  /* *** don't forget,  tail -> head */
  for(Edge i = 1; i <= nedges; i++) {
    Edge j = i/2;
    Vertex tail = tails[j];
    Vertex head = heads[j];
    tails[j] = tails[i-1];
    heads[j] = heads[i-1];
    tails[i-1] = tail;
    heads[i-1] = head;
  }
}

/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */

/*****************
 int SetEdge

 Set an edge value: set it to its new weight. Create if it
does not exist, destroy by setting to 0. 
*****************/
void SetEdge (Vertex tail, Vertex head, unsigned int weight, Network *nwp) 
{
  ENSURE_TH_ORDER;

  if(weight==0){
    DeleteEdgeFromTrees(tail,head,nwp);
  }else{
    if(EdgetreeSearch(tail, head, nwp->outedges)==0) AddEdgeToTrees(tail,head,nwp);
  }
}

Network *Redgelist2Network(SEXP elR, Rboolean empty){
  Vertex e = length(VECTOR_ELT(elR, 0));
  Vertex *tails = empty ? NULL : (Vertex*) INTEGER(VECTOR_ELT(elR, 0));
  Vertex *heads = empty ? NULL : (Vertex*) INTEGER(VECTOR_ELT(elR, 1));
  Vertex n = asInteger(getAttrib(elR, install("n")));
  Rboolean directed = asLogical(getAttrib(elR, install("directed")));
  Vertex bipartite = asInteger(getAttrib(elR, install("bipartite")));
  return NetworkInitialize(tails, heads, e, n, directed, bipartite, FALSE, 0, NULL); 
}


SEXP Network2Redgelist(Network *nwp){
  SEXP outl = PROTECT(mkNamed(VECSXP, (const char *[]){".tail", ".head", ""}));
  SEXP newnetworktails = PROTECT(allocVector(INTSXP, EDGECOUNT(nwp)));
  SEXP newnetworkheads = PROTECT(allocVector(INTSXP, EDGECOUNT(nwp)));
  EdgeTree2EdgeList((Vertex*)INTEGER(newnetworktails),
                    (Vertex*)INTEGER(newnetworkheads),
                    nwp,EDGECOUNT(nwp));
  SET_VECTOR_ELT(outl, 0, newnetworktails);
  SET_VECTOR_ELT(outl, 1, newnetworkheads);
  UNPROTECT(2);

  SEXP rownames = PROTECT(allocVector(INTSXP, EDGECOUNT(nwp)));
  int *r = INTEGER(rownames);
  for(unsigned int i=1; i<=EDGECOUNT(nwp); i++, r++) *r=i; 
  setAttrib(outl, install("row.names"), rownames);
  UNPROTECT(1);

  setAttrib(outl, install("n"), PROTECT(ScalarInteger(nwp->nnodes)));
  setAttrib(outl, install("directed"), PROTECT(ScalarLogical(nwp->directed_flag)));
  setAttrib(outl, install("bipartite"), PROTECT(ScalarInteger(nwp->bipartite)));
  UNPROTECT(3);

  SEXP class = PROTECT(mkRStrVec((const char*[]){"tibble_edgelist", "edgelist", "tbl_df", "tbl", "data.frame", NULL}));
  classgets(outl, class);
  UNPROTECT(1);

  UNPROTECT(1);
  return outl;
}
