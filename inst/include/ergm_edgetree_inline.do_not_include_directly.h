/*  File inst/include/ergm_edgetree_inline.do_not_include_directly.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2019 Statnet Commons
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
 Edge EdgetreePre(order)Successor

 Return the index of the next TreeNode in a preorder traversal.
*****************/
static inline Edge EdgetreePreSuccessor (TreeNode *edges, Edge x) {
  TreeNode *ptr;
  Edge y, z;

  // If we can go left, go left.
  if ((y=(ptr=edges+x)->left) != 0) 
    return y;
  // If we can go right, go right.
  if ((y=ptr->right) != 0) 
    return y;
  // Otherwise, keep going up until we can go right unless we just
  // went up from the right node.
  while ((y=ptr->parent)!=0){
    if((z=(ptr=(edges+y))->right)!=0 && z!=x) return z;
    x=y;
  }
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
  ENSURE_TH_ORDER;

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
 int GetEdge

Get edge value. Return 0 if edge does not exist.
*****************/
static inline unsigned int GetEdge (Vertex tail, Vertex head, Network *nwp) 
{
  ENSURE_TH_ORDER;

  return EdgetreeSearch(tail,head,nwp->outedges)!=0;
}

