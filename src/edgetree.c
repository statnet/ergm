/*  File src/edgetree.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2020 Statnet Commons
 */
#include "ergm_edgetree.h"

/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */


/*******************
 Network NetworkInitialize

 Initialize, construct binary tree version of network.  Note
 that the 0th TreeNode in the array is unused and should 
 have all its values set to zero
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

  if(lasttoggle_flag){
    nwp->duration_info.time=time;
    if(lasttoggle){
        nwp->duration_info.lasttoggle = (int *) Calloc(DYADCOUNT(nnodes, bipartite, directed_flag), int);
        memcpy(nwp->duration_info.lasttoggle, lasttoggle, DYADCOUNT(nnodes, bipartite, directed_flag) * sizeof(int));
    } else nwp->duration_info.lasttoggle = NULL;
  }
  else nwp->duration_info.lasttoggle = NULL;

  /*Configure a Network*/
  nwp->nnodes = nnodes;
  EDGECOUNT(nwp) = 0; /* Edges will be added one by one */
  nwp->directed_flag=directed_flag;
  nwp->bipartite=bipartite;

  DetShuffleEdges(tails,heads,nedges); /* shuffle to avoid worst-case performance */

  for(Edge i = 0; i < nedges; i++) {
    Vertex tail=tails[i], head=heads[i];
    if (!directed_flag && tail > head) 
      AddEdgeToTrees(head,tail,nwp); /* Undir edges always have tail < head */ 
    else 
      AddEdgeToTrees(tail,head,nwp);
  }
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
  Free(nwp->indegree);
  Free(nwp->outdegree);
  Free(nwp->inedges);
  Free(nwp->outedges);
  if(nwp->duration_info.lasttoggle){
    Free(nwp->duration_info.lasttoggle);
    nwp->duration_info.lasttoggle=NULL;
  }
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

  int directed_flag = dest->directed_flag = src->directed_flag;
  Vertex bipartite = dest->bipartite = src->bipartite;

  if(src->duration_info.lasttoggle){
    dest->duration_info.time=src->duration_info.time;
    dest->duration_info.lasttoggle = (int *) Calloc(DYADCOUNT(nnodes, bipartite, directed_flag), int);
    memcpy(dest->duration_info.lasttoggle, src->duration_info.lasttoggle,DYADCOUNT(nnodes, bipartite, directed_flag) * sizeof(int));
  }
  else dest->duration_info.lasttoggle = NULL;

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
  if (AddEdgeToTrees(tail,head,nwp))
    return 1;
  else 
    return 1 - DeleteEdgeFromTrees(tail,head,nwp);
}



/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */

/*****************
 Edge ToggleEdgeWithTimestamp
 By MSH 11/26/06

 Same as ToggleEdge, but this time with the additional
 step of updating the matrix of 'lasttoggle' times
 *****************/

/* *** don't forget tail->head, so this function now accepts tail before head */

int ToggleEdgeWithTimestamp(Vertex tail, Vertex head, Network *nwp){
  Edge k;

  /* don't forget, tails < heads in undirected networks now  */
  ENSURE_TH_ORDER;
  
  if(nwp->duration_info.lasttoggle){ /* Skip timestamps if no duration info. */
    if(nwp->bipartite){
      k = (head-nwp->bipartite-1)*(nwp->bipartite) + tail - 1;
    }else{
      if (nwp->directed_flag) 
	k = (head-1)*(nwp->nnodes-1) + tail - ((tail>head) ? 1:0) - 1; 
      else
	k = (head-1)*(head-2)/2 + tail - 1;    
    }
    nwp->duration_info.lasttoggle[k] = nwp->duration_info.time;
  }
  
  if (AddEdgeToTrees(tail,head,nwp))
    return 1;
  else 
    return 1 - DeleteEdgeFromTrees(tail,head,nwp);
}

/*****************
 void TouchEdge

 Named after the UNIX "touch" command.
 Set an edge's time-stamp to the current MCMC time.
 *****************/

/* *** don't forget tail->head, so this function now accepts tail before head */

void TouchEdge(Vertex tail, Vertex head, Network *nwp){
  unsigned int k;
  if(nwp->duration_info.lasttoggle){ /* Skip timestamps if no duration info. */
    if(nwp->bipartite){
      k = (head-nwp->bipartite-1)*(nwp->bipartite) + tail - 1;
    }else{
      if (nwp->directed_flag) 
	k = (head-1)*(nwp->nnodes-1) + tail - ((tail>head) ? 1:0) - 1; 
      else
	k = (head-1)*(head-2)/2 + tail - 1;    
    }
    nwp->duration_info.lasttoggle[k] = nwp->duration_info.time;
  }
}


/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */

/*****************
 Edge AddEdgeToTrees

 Add an edge from tail to head after checking to see
 if it's legal. Return 1 if edge added, 0 otherwise.  Since each
 "edge" should be added to both the list of outedges and the list of 
 inedges, this actually involves two calls to AddHalfedgeToTree (hence
 "Trees" instead of "Tree" in the name of this function).
*****************/

/* *** don't forget tail->head, so this function now accepts tail before head */

int AddEdgeToTrees(Vertex tail, Vertex head, Network *nwp){
  if (EdgetreeSearch(tail, head, nwp->outedges) == 0) {
    AddHalfedgeToTree(tail, head, nwp->outedges, &(nwp->last_outedge));
    AddHalfedgeToTree(head, tail, nwp->inedges, &(nwp->last_inedge));
    ++nwp->outdegree[tail];
    ++nwp->indegree[head];
    ++EDGECOUNT(nwp);
    CheckEdgetreeFull(nwp); 
    return 1;
  }
  return 0;
}

/*****************
 void AddHalfedgeToTree:  Only called by AddEdgeToTrees
*****************/
void AddHalfedgeToTree (Vertex a, Vertex b, TreeNode *edges, Edge *last_edge){
  TreeNode *eptr = edges+a, *newnode;
  Edge e;

  if (eptr->value==0) { /* This is the first edge for vertex a. */
    eptr->value=b;
    return;
  }
  (newnode = edges + (++*last_edge))->value=b;  
  newnode->left = newnode->right = 0;
  /* Now find the parent of this new edge */
  for (e=a; e!=0; e=(b < (eptr=edges+e)->value) ? eptr->left : eptr->right);
  newnode->parent=eptr-edges;  /* Point from the new edge to the parent... */
  if (b < eptr->value)  /* ...and have the parent point back. */
    eptr->left=*last_edge; 
  else
    eptr->right=*last_edge;
}

/*****************
void CheckEdgetreeFull
*****************/
void CheckEdgetreeFull (Network *nwp) {
  const unsigned int mult=2;
  
  // Note that maximum index in the nwp->*edges is nwp->maxedges-1, and we need to keep one element open for the next insertion.
  if(nwp->last_outedge==nwp->maxedges-2 || nwp->last_inedge==nwp->maxedges-2){
    // Only enlarge the non-root part of the array.
    Edge newmax = nwp->nnodes + 1 + (nwp->maxedges - nwp->nnodes - 1)*mult;
    nwp->inedges = (TreeNode *) Realloc(nwp->inedges, newmax, TreeNode);
    memset(nwp->inedges+nwp->maxedges, 0,
	   sizeof(TreeNode) * (newmax-nwp->maxedges));
    nwp->outedges = (TreeNode *) Realloc(nwp->outedges, newmax, TreeNode);
    memset(nwp->outedges+nwp->maxedges, 0,
	   sizeof(TreeNode) * (newmax-nwp->maxedges));
    nwp->maxedges = newmax;
  }
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
  if (DeleteHalfedgeFromTree(tail, head, nwp->outedges,&(nwp->last_outedge))&&
      DeleteHalfedgeFromTree(head, tail, nwp->inedges, &(nwp->last_inedge))) {
    --nwp->outdegree[tail];
    --nwp->indegree[head];
    --EDGECOUNT(nwp);
    if(nwp->last_outedge < nwp->nnodes) nwp->last_outedge=nwp->nnodes;
    if(nwp->last_inedge < nwp->nnodes) nwp->last_inedge=nwp->nnodes;
    return 1;
  }
  return 0;
}

/*****************
 int DeleteHalfedgeFromTree

 Delete the TreeNode with value b from the tree rooted at edges[a].
 Return 0 if no such TreeNode exists, 1 otherwise.  Also update the
 value of *last_edge appropriately.
*****************/
int DeleteHalfedgeFromTree(Vertex a, Vertex b, TreeNode *edges,
		     Edge *last_edge){
  Edge x, z, root=(Edge)a;
  TreeNode *xptr, *zptr, *ptr;

  if ((z=EdgetreeSearch(a, b, edges))==0)  /* z is the current TreeNode. */
    return 0;  /* This edge doesn't exist, so return 0 */
  /* First, determine which node to splice out; this is z.  If the current
     z has two children, then we'll actually splice out its successor. */
  if ((zptr=edges+z)->left != 0 && zptr->right != 0) {
    /* Select which child to promote based on whether the left child's
       position is divisible by 2: the position of a node in an edge
       tree is effectively random, *unless* it's a root node. Using
       the left child ensures that it is not a root node. */
    if(zptr->left&1u)
      z=EdgetreeSuccessor(edges, z);  
    else
      z=EdgetreePredecessor(edges, z);  
    zptr->value = (ptr=edges+z)->value;
    zptr=ptr;
  }
  /* Set x to the child of z (there is at most one). */
  if ((x=zptr->left) == 0)
    x = zptr->right;
  /* Splice out node z */
  if (z == root) {
    zptr->value = (xptr=edges+x)->value;
    if (x != 0) {
      if ((zptr->left=xptr->left) != 0)
	(edges+zptr->left)->parent = z;
      if ((zptr->right=xptr->right) != 0)
	(edges+zptr->right)->parent = z;
      zptr=edges+(z=x);
    }  else 
      return 1;
  } else {
    if (x != 0)
      (xptr=edges+x)->parent = zptr->parent;
    if (z==(ptr=(edges+zptr->parent))->left)
      ptr->left = x;
    else 
      ptr->right = x;
  }  
  /* Clear z node, update *last_edge if necessary. */
  zptr->value=0;
  if(z!=root){
    RelocateHalfedge(*last_edge,z,edges);
    (*last_edge)--;
  }
  return 1;
}

void RelocateHalfedge(Edge from, Edge to, TreeNode *edges){
  if(from==to) return;
  TreeNode *toptr=edges+to, *fromptr=edges+from;

  if(fromptr->left) edges[fromptr->left].parent = to;
  if(fromptr->right) edges[fromptr->right].parent = to;
  if(fromptr->parent){
    TreeNode *parentptr = edges+fromptr->parent;
    if(parentptr->left==from) parentptr->left = to;
    else parentptr->right =  to;
  }
  memcpy(toptr,fromptr,sizeof(TreeNode));
  fromptr->value = 0;
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
    AddEdgeToTrees(tail,head,nwp);
  }
}

/*****************
 Edge SetEdgeWithTimestamp

 Same as SetEdge, but this time with the additional
 step of updating the matrix of 'lasttoggle' times
 *****************/
void SetEdgeWithTimestamp (Vertex tail, Vertex head, unsigned int weight, Network *nwp) 
{
  Edge k;

  ENSURE_TH_ORDER;
  
  if(nwp->duration_info.lasttoggle){ /* Skip timestamps if no duration info. */
    if(nwp->bipartite){
      k = (head-nwp->bipartite-1)*(nwp->bipartite) + tail - 1;
    }else{
      if (nwp->directed_flag) 
	k = (head-1)*(nwp->nnodes-1) + tail - ((tail>head) ? 1:0) - 1; 
      else
	k = (head-1)*(head-2)/2 + tail - 1;    
    }
    nwp->duration_info.lasttoggle[k] = nwp->duration_info.time;
  }

  SetEdge(tail,head,weight,nwp);
}
