/*  File src/wtedgetree.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2018 Statnet Commons
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
			      int lasttoggle_flag, int time, int *lasttoggle) {
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
  PutRNGstate();
  return nw;
}


/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */

/*Takes vectors of doubles for edges; used only when constructing from inputparams. */
WtNetwork WtNetworkInitializeD(double *tails, double *heads, double *weights, Edge nedges,
           Vertex nnodes, int directed_flag, Vertex bipartite,
           int lasttoggle_flag, int time, int *lasttoggle) {

  /* *** don't forget, tail -> head */

  Vertex *itails=(Vertex*)malloc(sizeof(Vertex)*nedges);
  Vertex *iheads=(Vertex*)malloc(sizeof(Vertex)*nedges);
  
  for(Edge i=0; i<nedges; i++){
    itails[i]=tails[i];
    iheads[i]=heads[i];
  }

  WtNetwork nw=WtNetworkInitialize(itails,iheads,weights,nedges,nnodes,directed_flag,bipartite,lasttoggle_flag, time, lasttoggle);

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

  return dest;
}

/*****************
 Edge WtEdgetreeSearch

 Check to see if there's a WtTreeNode with value b 
 in the tree rooted at edges[a].  Return i such that 
 edges[i] is that WtTreeNode, or 0 if none.
*****************/
Edge WtEdgetreeSearch (Vertex a, Vertex b, WtTreeNode *edges) {
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
 Edge EdgetreeSuccessor

 Return the index of the WtTreeNode with the smallest value
 greater than edges[x].value in the same edge tree, or 0
 if none.  This is used by (for instance)
 the DeleteHalfedgeFromTree function.
*****************/
Edge WtEdgetreeSuccessor (WtTreeNode *edges, Edge x) {
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
Edge WtEdgetreePredecessor (WtTreeNode *edges, Edge x) {
  WtTreeNode *ptr;
  Edge y;

  if ((y=(ptr=edges+x)->left) != 0) 
    return WtEdgetreeMaximum (edges, y);
  while ((y=ptr->parent)!=0 && x==(ptr=edges+y)->left) 
    x=y;
  return y; 
}   

/*****************
 Edge WtEdgetreeMinimum

 Return the index of the WtTreeNode with the
 smallest value in the subtree rooted at x
*****************/
Edge WtEdgetreeMinimum (WtTreeNode *edges, Edge x) {
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
Edge WtEdgetreeMaximum (WtTreeNode *edges, Edge x) {
  Edge y;

  while ((y=(edges+x)->right) != 0)
    x=y;
  return x;
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

int WtToggleEdgeWithTimestamp (Vertex tail, Vertex head, double weight, WtNetwork *nwp) 
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
void WtSetEdge (Vertex tail, Vertex head, double weight, WtNetwork *nwp) 
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
double WtGetEdge (Vertex tail, Vertex head, WtNetwork *nwp) 
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
void WtSetEdgeWithTimestamp (Vertex tail, Vertex head, double weight, WtNetwork *nwp) 
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

/*****************
 long int ElapsedTime

 Return time since given (tail,head) was last toggled using
 ToggleEdgeWithTimestamp function
*****************/

/* *** don't forget tail->head, so this function now accepts tail before head */

int WtElapsedTime (Vertex tail, Vertex head, WtNetwork *nwp){
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



/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */


/*****************
 void TouchEdge

 Named after the UNIX "touch" command.
 Set an edge's time-stamp to the current MCMC time.
 *****************/

/* *** don't forget tail->head, so this function now accepts tail before head */

void WtTouchEdge(Vertex tail, Vertex head, WtNetwork *nwp){
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

int WtAddEdgeToTrees(Vertex tail, Vertex head, double weight, WtNetwork *nwp){
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

/*****************
 void WtAddHalfedgeToTree:  Only called by WtAddEdgeToTrees
*****************/
void WtAddHalfedgeToTree (Vertex a, Vertex b, double weight, WtTreeNode *edges, Edge *last_edge){
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
void CheckEdgetreeFull
*****************/
void WtCheckEdgetreeFull (WtNetwork *nwp) {
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
 int WtDeleteHalfedgeFromTree

 Delete the WtTreeNode with value b from the tree rooted at edges[a].
 Return 0 if no such WtTreeNode exists, 1 otherwise.  Also update the
 value of *last_edge appropriately.
*****************/
int WtDeleteHalfedgeFromTree(Vertex a, Vertex b, WtTreeNode *edges,
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

void WtRelocateHalfedge(Edge from, Edge to, WtTreeNode *edges){
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
  
  // If the index is too high or too low, exit immediately.
  if (i > ndyads - nwp->nedges || i<=0)
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
