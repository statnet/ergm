/*  File src/wtedgetree_inline.do_not_include_directly.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2021 Statnet Commons
 */

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

/*****************
 int WtDeleteHalfedgeFromTree

 Delete the WtTreeNode with value b from the tree rooted at edges[a].
 Return 0 if no such WtTreeNode exists, 1 otherwise.  Also update the
 value of *last_edge appropriately.
*****************/
static inline void WtDeleteHalfedgeFromTreeAt(Vertex a, Vertex b, WtTreeNode *edges,
                                              Edge *last_edge, Edge z){ 
  Edge x, root=(Edge)a;
  WtTreeNode *xptr, *zptr, *ptr;

  /* First, determine which node to splice out; this is z.  If the current
     z has two children, then we'll actually splice out its successor. */
  if ((zptr=edges+z)->left != 0 && zptr->right != 0) {
    /* Select which child to promote based on whether the left child's
       position is divisible by 2: the position of a node in an edge
       tree is effectively random, *unless* it's a root node. Using
       the left child ensures that it is not a root node. */
    if(zptr->left&1u)
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
      return;
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
  return;
}

/*****************
void CheckEdgetreeFull
*****************/
static inline void WtCheckEdgetreeFull (WtNetwork *nwp) {
  const unsigned int mult=2;
  
  // Note that maximum index in the nwp->*edges is nwp->maxedges-1, and we need to keep one element open for the next insertion.
  if(nwp->last_outedge==nwp->maxedges-2 || nwp->last_inedge==nwp->maxedges-2){
    // Only enlarge the non-root part of the array.
    Edge newmax = nwp->nnodes + 1 + (nwp->maxedges - nwp->nnodes - 1)*mult;
    nwp->inedges = (WtTreeNode *) Realloc(nwp->inedges, newmax, WtTreeNode);
    memset(nwp->inedges+nwp->maxedges, 0,
	   sizeof(WtTreeNode) * (newmax-nwp->maxedges));
    nwp->outedges = (WtTreeNode *) Realloc(nwp->outedges, newmax, WtTreeNode);
    memset(nwp->outedges+nwp->maxedges, 0,
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
