/*  File src/edgetree_inline.h in package ergm, part of the Statnet suite
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
 Edge EdgetreeSearch

 Check to see if there's a TreeNode with value b 
 in the tree rooted at edges[a].  Return i such that 
 edges[i] is that TreeNode, or 0 if none.
*****************/
static inline Edge EdgetreeSearch (Vertex a, Vertex b, TreeNode *edges) {
  TreeNode *es;
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
 Edge EdgetreeMinimum

 Return the index of the TreeNode with the
 smallest value in the subtree rooted at x
*****************/
static inline Edge EdgetreeMinimum (TreeNode *edges, Edge x) {
  Edge y;

  while ((y=(edges+x)->left) != 0)
    x=y;
  return x;
}

/*****************
 Edge EdgetreeMaximum

 Return the index of the TreeNode with the
 greatest value in the subtree rooted at x
*****************/
static inline Edge EdgetreeMaximum (TreeNode *edges, Edge x) {
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

 Return the index of the TreeNode with the smallest value
 greater than edges[x].value in the same edge tree, or 0
 if none.  This is used by (for instance)
 the DeleteHalfedgeFromTree function.
*****************/
static inline Edge EdgetreeSuccessor (TreeNode *edges, Edge x) {
  TreeNode *ptr;
  Edge y;

  if ((y=(ptr=edges+x)->right) != 0) 
    return EdgetreeMinimum (edges, y);
  while ((y=ptr->parent)!=0 && x==(ptr=edges+y)->right) 
    x=y;
  return y; 
}   

/*****************
 Edge EdgetreePredecessor

 Return the index of the TreeNode with the smallest value
 greater than edges[x].value in the same edge tree, or 0
 if none.  This is used by (for instance)
 the DeleteHalfedgeFromTree function.
*****************/
static inline Edge EdgetreePredecessor (TreeNode *edges, Edge x) {
  TreeNode *ptr;
  Edge y;

  if ((y=(ptr=edges+x)->left) != 0) 
    return EdgetreeMaximum (edges, y);
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

static inline int ElapsedTime(Vertex tail, Vertex head, Network *nwp){
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

static inline void TouchEdge(Vertex tail, Vertex head, Network *nwp){
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
static inline void CheckEdgetreeFull
*****************/
static inline void CheckEdgetreeFull (Network *nwp) {
  const unsigned int mult=2;
  
  // Note that maximum index in the nwp->*edges is nwp->maxedges-1, and we need to keep one element open for the next insertion.
  if(nwp->last_outedge==nwp->maxedges-2 || nwp->last_inedge==nwp->maxedges-2){
    // Only enlarge the non-root part of the array.
    Edge newmax = nwp->maxedges + (nwp->maxedges - nwp->nnodes - 1)*mult;
    nwp->inedges = (TreeNode *) realloc(nwp->inedges, 
					  sizeof(TreeNode) * newmax);
    memset(nwp->inedges+nwp->last_inedge+2,0,
	   sizeof(TreeNode) * (newmax-nwp->maxedges));
    nwp->outedges = (TreeNode *) realloc(nwp->outedges, 
					   sizeof(TreeNode) * newmax);
    memset(nwp->outedges+nwp->last_outedge+2,0,
	   sizeof(TreeNode) * (newmax-nwp->maxedges));
    nwp->maxedges = newmax;
  }
}

/*****************
 void AddHalfedgeToTree:  Only called by AddEdgeToTrees
*****************/
static inline void AddHalfedgeToTree (Vertex a, Vertex b, TreeNode *edges, Edge *last_edge){
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
 Edge AddEdgeToTrees

 Add an edge from tail to head after checking to see
 if it's legal. Return 1 if edge added, 0 otherwise.  Since each
 "edge" should be added to both the list of outedges and the list of 
 inedges, this actually involves two calls to AddHalfedgeToTree (hence
 "Trees" instead of "Tree" in the name of this function).
*****************/

/* *** don't forget tail->head, so this function now accepts tail before head */

static inline int AddEdgeToTrees(Vertex tail, Vertex head, Network *nwp){
  if (EdgetreeSearch(tail, head, nwp->outedges) == 0) {
    AddHalfedgeToTree(tail, head, nwp->outedges, &(nwp->last_outedge));
    AddHalfedgeToTree(head, tail, nwp->inedges, &(nwp->last_inedge));
    ++nwp->outdegree[tail];
    ++nwp->indegree[head];
    ++nwp->nedges;
    CheckEdgetreeFull(nwp); 
    return 1;
  }
  return 0;
}

static inline void RelocateHalfedge(Edge from, Edge to, TreeNode *edges){
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


/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */

/*****************
 int DeleteHalfedgeFromTree

 Delete the TreeNode with value b from the tree rooted at edges[a].
 Return 0 if no such TreeNode exists, 1 otherwise.  Also update the
 value of *last_edge appropriately.
*****************/
static inline int DeleteHalfedgeFromTree(Vertex a, Vertex b, TreeNode *edges,
		     Edge *last_edge){
  Edge x, z, root=(Edge)a;
  TreeNode *xptr, *zptr, *ptr;

  if ((z=EdgetreeSearch(a, b, edges))==0)  /* z is the current TreeNode. */
    return 0;  /* This edge doesn't exist, so return 0 */
  /* First, determine which node to splice out; this is z.  If the current
     z has two children, then we'll actually splice out its successor. */
  if ((zptr=edges+z)->left != 0 && zptr->right != 0) {
    if(unif_rand()<0.5)
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

/*****************
 int DeleteEdgeFromTrees

 Find and delete the edge from tail to head.  
 Return 1 if successful, 0 otherwise.  As with AddEdgeToTrees, this must
 be done once for outedges and once for inedges.
*****************/

/* *** don't forget tail->head, so this function now accepts tail before head */

static inline int DeleteEdgeFromTrees(Vertex tail, Vertex head, Network *nwp){
  if (DeleteHalfedgeFromTree(tail, head, nwp->outedges,&(nwp->last_outedge))&&
      DeleteHalfedgeFromTree(head, tail, nwp->inedges, &(nwp->last_inedge))) {
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
 Edge ToggleEdge

 Toggle an edge:  Set it to the opposite of its current
 value.  Return 1 if edge added, 0 if deleted.
*****************/

/* *** don't forget tail->head, so this function now accepts tail before head */

static inline int ToggleEdge (Vertex tail, Vertex head, Network *nwp) 
{
  /* don't forget tails < heads now for undirected networks */
  if (!(nwp->directed_flag) && tail > head) {
    Vertex temp;
    temp = tail; /*  Make sure tail<head always for undirected edges */
    tail = head;
    head = temp;
  }
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

static inline int ToggleEdgeWithTimestamp(Vertex tail, Vertex head, Network *nwp){
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
  
  if (AddEdgeToTrees(tail,head,nwp))
    return 1;
  else 
    return 1 - DeleteEdgeFromTrees(tail,head,nwp);
}

/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */



