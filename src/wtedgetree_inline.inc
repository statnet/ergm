/*  File src/wtedgetree_inline.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2013 Statnet Commons
 */

/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */

/*****************
 Edge WtEdgetreeSearch

 Check to see if there's a WtTreeNode with value b 
 in the tree rooted at edges[a].  Return i such that 
 edges[i] is that WtTreeNode, or 0 if none.
*****************/
static inline Edge WtEdgetreeSearch (Vertex a, Vertex b, WtTreeNode *edges) {
  WtTreeNode *es;
  Edge e = a;
  Vertex v;

  es = edges + e;
  v = es->value;
  while(e != 0 && b != v)  {
      e = (b<v)?  es->left : es->right;
      es = edges + e;
      v = es->value;
    }
  return e;
}

/*****************
 Edge WtEdgetreeMinimum

 Return the index of the WtTreeNode with the
 smallest value in the subtree rooted at x
*****************/
static inline Edge WtEdgetreeMinimum (WtTreeNode *edges, Edge x) {
  Edge y;

  while ((y=(edges+x)->left) != 0)
    x=y;
  return x;
}

/*****************
 Edge WtEdgetreeMaximum

 Return the index of the WtTreeNode with the
 greatest value in the subtree rooted at x
*****************/
static inline Edge WtEdgetreeMaximum (WtTreeNode *edges, Edge x) {
  Edge y;

  while ((y=(edges+x)->right) != 0)
    x=y;
  return x;
}

/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */

/*****************
 Edge EdgetreeSuccessor

 Return the index of the WtTreeNode with the smallest value
 greater than edges[x].value in the same edge tree, or 0
 if none.  This is used by (for instance)
 the DeleteHalfedgeFromTree function.
*****************/
static inline Edge WtEdgetreeSuccessor (WtTreeNode *edges, Edge x) {
  WtTreeNode *ptr;
  Edge y;

  if ((y=(ptr=edges+x)->right) != 0) 
    return WtEdgetreeMinimum (edges, y);
  while ((y=ptr->parent)!=0 && x==(ptr=edges+y)->right) 
    x=y;
  return y; 
}   

/*****************
 Edge WtEdgetreePredecessor

 Return the index of the WtTreeNode with the smallest value
 greater than edges[x].value in the same edge tree, or 0
 if none.  This is used by (for instance)
 the WtDeleteHalfedgeFromTree function.
*****************/
static inline Edge WtEdgetreePredecessor (WtTreeNode *edges, Edge x) {
  WtTreeNode *ptr;
  Edge y;

  if ((y=(ptr=edges+x)->left) != 0) 
    return WtEdgetreeMaximum (edges, y);
  while ((y=ptr->parent)!=0 && x==(ptr=edges+y)->left) 
    x=y;
  return y; 
}   

/*****************
 long int ElapsedTime

 Return time since given (tail,head) was last toggled using
 ToggleEdgeWithTimestamp function
*****************/

/* *** don't forget tail->head, so this function now accepts tail before head */

static inline int WtElapsedTime (Vertex tail, Vertex head, WtNetwork *nwp){
  Edge k;
  /* don't forget, tails < heads now in undirected networks */
  if (!(nwp->directed_flag) && tail > head) {
    Vertex temp;
    temp = tail; /*  Make sure tail<head always for undirected edges */
    tail = head;
    head = temp;
  }

  if(nwp->duration_info.lasttoggle){ /* Return 0 if no duration info. */
    if(nwp->bipartite){
      k = (head-nwp->bipartite-1)*(nwp->bipartite) + tail - 1;
    }else{
      if (nwp->directed_flag) 
	k = (head-1)*(nwp->nnodes-1) + tail - ((tail>head) ? 1:0) - 1; 
      else
	k = (head-1)*(head-2)/2 + tail - 1;    
    }
    return nwp->duration_info.time - nwp->duration_info.lasttoggle[k];
  }
  else return 0; 
  /* Should maybe return an error code of some sort, since 0 elapsed time
     is valid output. Need to think about it. */
}


/*****************
 void TouchEdge

 Named after the UNIX "touch" command.
 Set an edge's time-stamp to the current MCMC time.
 *****************/

/* *** don't forget tail->head, so this function now accepts tail before head */

static inline void WtTouchEdge(Vertex tail, Vertex head, WtNetwork *nwp){
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
void CheckEdgetreeFull
*****************/
static inline void WtCheckEdgetreeFull (WtNetwork *nwp) {
  const unsigned int mult=2;
  
  // Note that maximum index in the nwp->*edges is nwp->maxedges-1, and we need to keep one element open for the next insertion.
  if(nwp->last_outedge==nwp->maxedges-2 || nwp->last_inedge==nwp->maxedges-2){
    // Only enlarge the non-root part of the array.
    Edge newmax = nwp->maxedges + (nwp->maxedges - nwp->nnodes - 1)*mult;
    nwp->inedges = (WtTreeNode *) realloc(nwp->inedges, 
					  sizeof(WtTreeNode) * newmax);
    memset(nwp->inedges+nwp->last_inedge+2,0,
	   sizeof(WtTreeNode) * (newmax-nwp->maxedges));
    nwp->outedges = (WtTreeNode *) realloc(nwp->outedges, 
					   sizeof(WtTreeNode) * newmax);
    memset(nwp->outedges+nwp->last_outedge+2,0,
	   sizeof(WtTreeNode) * (newmax-nwp->maxedges));
    nwp->maxedges = newmax;
  }
}

/*****************
 void WtAddHalfedgeToTree:  Only called by WtAddEdgeToTrees
*****************/
static inline void WtAddHalfedgeToTree (Vertex a, Vertex b, double weight, WtTreeNode *edges, Edge *last_edge){
  WtTreeNode *eptr = edges+a, *newnode;
  Edge e;

  if (eptr->value==0) { /* This is the first edge for vertex a. */
    eptr->value=b;
    eptr->weight = weight;  /*  Add weight too */
    return;
  }
  (newnode = edges + (++*last_edge))->value=b;  
  newnode->left = newnode->right = 0;
  newnode->weight=weight;  /*  Add weight too */
  /* Now find the parent of this new edge */
  for (e=a; e!=0; e=(b < (eptr=edges+e)->value) ? eptr->left : eptr->right);
  newnode->parent=eptr-edges;  /* Point from the new edge to the parent... */
  if (b < eptr->value)  /* ...and have the parent point back. */
    eptr->left=*last_edge; 
  else
    eptr->right=*last_edge;
}


/*****************
 Edge AddEdgeToTrees

 Add an edge from tail to head after checking to see
 if it's legal. Return 1 if edge added, 0 otherwise.  Since each
 "edge" should be added to both the list of outedges and the list of 
 inedges, this actually involves two calls to AddHalfedgeToTree (hence
 "Trees" instead of "Tree" in the name of this function).
*****************/

/* *** don't forget tail->head, so this function now accepts tail before head */

static inline int WtAddEdgeToTrees(Vertex tail, Vertex head, double weight, WtNetwork *nwp){
  if (WtEdgetreeSearch(tail, head, nwp->outedges) == 0) {
    WtAddHalfedgeToTree(tail, head, weight, nwp->outedges, &(nwp->last_outedge));
    WtAddHalfedgeToTree(head, tail, weight, nwp->inedges, &(nwp->last_inedge));
    ++nwp->outdegree[tail];
    ++nwp->indegree[head];
    ++nwp->nedges;
    WtCheckEdgetreeFull(nwp); 
    return 1;
  }
  return 0;
}

static inline void WtRelocateHalfedge(Edge from, Edge to, WtTreeNode *edges){
  if(from==to) return;
  WtTreeNode *toptr=edges+to, *fromptr=edges+from;

  if(fromptr->left) edges[fromptr->left].parent = to;
  if(fromptr->right) edges[fromptr->right].parent = to;
  if(fromptr->parent){
    WtTreeNode *parentptr = edges+fromptr->parent;
    if(parentptr->left==from) parentptr->left = to;
    else parentptr->right =  to;
  }
  memcpy(toptr,fromptr,sizeof(WtTreeNode));
  fromptr->value = 0;
}


/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */

/*****************
 int WtDeleteHalfedgeFromTree

 Delete the WtTreeNode with value b from the tree rooted at edges[a].
 Return 0 if no such WtTreeNode exists, 1 otherwise.  Also update the
 value of *last_edge appropriately.
*****************/
static inline int WtDeleteHalfedgeFromTree(Vertex a, Vertex b, WtTreeNode *edges,
		     Edge *last_edge){ 
  Edge x, z, root=(Edge)a;
  WtTreeNode *xptr, *zptr, *ptr;

  if ((z=WtEdgetreeSearch(a, b, edges))==0)  /* z is the current WtTreeNode. */
    return 0; /* This edge doesn't exist, so return 0 */
  /* First, determine which node to splice out; this is z.  If the current
     z has two children, then we'll actually splice out its successor. */
  if ((zptr=edges+z)->left != 0 && zptr->right != 0) {
    if(unif_rand()<0.5)
      z=WtEdgetreeSuccessor(edges, z);  
    else
      z=WtEdgetreePredecessor(edges, z);  
    zptr->value = (ptr=edges+z)->value;
    zptr->weight = ptr->weight;
    zptr=ptr;
  }
  /* Set x to the child of z (there is at most one). */
  if ((x=zptr->left) == 0)
    x = zptr->right;
  /* Splice out node z */
  if (z == root) {
    zptr->value = (xptr=edges+x)->value;
    zptr->weight = xptr->weight;
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
    WtRelocateHalfedge(*last_edge,z,edges);
    (*last_edge)--;
  }
  return 1;
}

/*****************
 int WtDeleteEdgeFromTrees

 Find and delete the edge from tail to head.  
 Return 1 if successful, 0 otherwise.  As with AddEdgeToTrees, this must
 be done once for outedges and once for inedges.
*****************/

/* *** don't forget tail->head, so this function now accepts tail before head */

static inline int WtDeleteEdgeFromTrees(Vertex tail, Vertex head, WtNetwork *nwp){
  if (WtDeleteHalfedgeFromTree(tail, head, nwp->outedges,&(nwp->last_outedge))&&
      WtDeleteHalfedgeFromTree(head, tail, nwp->inedges, &(nwp->last_inedge))) {
    --nwp->outdegree[tail];
    --nwp->indegree[head];
    --nwp->nedges;
    if(nwp->last_outedge < nwp->nnodes) nwp->last_outedge=nwp->nnodes;
    if(nwp->last_inedge < nwp->nnodes) nwp->last_inedge=nwp->nnodes;
    return 1;
  }
  return 0;
}


/*****************
 int WtToggleEdge

 Toggle an edge:  Set it to the opposite of its current
 value.  Return 1 if edge added, 0 if deleted.
*****************/

/* *** don't forget tail->head, so this function now accepts tail before head */

static inline int WtToggleEdge (Vertex tail, Vertex head, double weight, WtNetwork *nwp) 
{
  /* don't forget tails < heads now for undirected networks */
  if (!(nwp->directed_flag) && tail > head) {
    Vertex temp;
    temp = tail; /*  Make sure tail<head always for undirected edges */
    tail = head;
    head = temp;
  }
  if (WtAddEdgeToTrees(tail,head,weight,nwp))
    return 1;
  else 
    return 1 - WtDeleteEdgeFromTrees(tail,head,nwp);
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

static inline int WtToggleEdgeWithTimestamp (Vertex tail, Vertex head, double weight, WtNetwork *nwp) 
{
  Edge k;

  /* don't forget, tails < heads in undirected networks now  */
  if (!(nwp->directed_flag) && tail > head) {
    Vertex temp;
    temp = tail; /*  Make sure tail<head always for undirected edges */
    tail = head;
    head = temp;
  }
  
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

  if (WtAddEdgeToTrees(tail,head,weight,nwp))
    return 1;
  else 
    return 1 - WtDeleteEdgeFromTrees(tail,head,nwp);
}

/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */

/*****************
 int WtSetEdge

 Set an weighted edge value: set it to its new weight. Create if it
does not exist, destroy by setting to 0. 
*****************/
static inline void WtSetEdge (Vertex tail, Vertex head, double weight, WtNetwork *nwp) 
{
  if (!(nwp->directed_flag) && tail > head) {
    Vertex temp;
    temp = tail; /*  Make sure tail<head always for undirected edges */
    tail = head;
    head = temp;
  }

  if(weight==0){
    // If the function is to set the edge value to 0, just delete it.
    WtDeleteEdgeFromTrees(tail,head,nwp);
  }else{
    // Find the out-edge
    Edge oe=WtEdgetreeSearch(tail,head,nwp->outedges);
    if(oe){
      // If it exists AND already has the target weight, do nothing.
      if(nwp->outedges[oe].weight==weight) return;
      else{
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

/*****************
 int WtGetEdge

Get weighted edge value. Return 0 if edge does not exist.
*****************/
static inline double WtGetEdge (Vertex tail, Vertex head, WtNetwork *nwp) 
{
  if (!(nwp->directed_flag) && tail > head) {
    Vertex temp;
    temp = tail; /*  Make sure tail<head always for undirected edges */
    tail = head;
    head = temp;
  }

  Edge oe=WtEdgetreeSearch(tail,head,nwp->outedges);
  if(oe) return nwp->outedges[oe].weight;
  else return 0;
}

/*****************
 Edge WtSetEdgeWithTimestamp

 Same as WtSetEdge, but this time with the additional
 step of updating the matrix of 'lasttoggle' times
 *****************/
static inline void WtSetEdgeWithTimestamp (Vertex tail, Vertex head, double weight, WtNetwork *nwp) 
{
  Edge k;

  if (!(nwp->directed_flag) && tail > head) {
    Vertex temp;
    temp = tail; /*  Make sure tail<head always for undirected edges */
    tail = head;
    head = temp;
  }
  
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

  WtSetEdge(tail,head,weight,nwp);
}



/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */

