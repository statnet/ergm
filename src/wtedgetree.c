#include "wtedgetree.h"

/*######################################################################
#
# copyright (c) 2003, Mark S. Handcock, University of Washington
#                     David R. Hunter, Penn State University
#                     Carter T. Butts, University of California - Irvine
#                     Martina Morris, University of Washington
# 
# For license information see http://statnetproject.org/license
#
# We have invested a lot of time and effort in creating 'statnet',
# for use by other researchers. We require that the attributions
# in the software are retained (even if only pieces of it are used),
# and that there is attribution when the package is loaded (e.g., via
# "library" or "require"). This is to stop
# "rebadging" of the software. 
#
# Cite us!
#
# To cite see http://statnetproject.org/cite
#
# .First.lib is run when the package is loaded.
#
###################################################################### */

/*******************
 WtNetwork WtNetworkInitialize

 Initialize, construct binary tree version of network with weights.  Note
 that the 0th WtTreeNode in the array is unused and should 
 have all its values set to zero
*******************/

/* *** don't forget tail->head, so this function now accepts tails before heads */

WtNetwork WtNetworkInitialize(Vertex *tails, Vertex *heads, double *weights,
			      Edge nedges, Vertex nnodes, int directed_flag, Vertex bipartite,
			      int lasttoggle_flag) {
  Edge i;
  WtNetwork nw;

  GetRNGstate();  /* R function enabling uniform RNG */

  nw.next_inedge = nw.next_outedge = (Edge)nnodes+1;
  /* Calloc will zero the allocated memory for us, probably a lot
     faster. */
  nw.outdegree = (Vertex *) calloc((nnodes+1),sizeof(Vertex));
  nw.indegree  = (Vertex *) calloc((nnodes+1),sizeof(Vertex));
  nw.maxedges = MAX(nedges,1)+nnodes+2; /* Maybe larger than needed? */
  nw.inedges = (WtTreeNode *) calloc(nw.maxedges,sizeof(WtTreeNode));
  nw.outedges = (WtTreeNode *) calloc(nw.maxedges,sizeof(WtTreeNode));

  if(lasttoggle_flag){
    nw.duration_info.MCMCtimer=0;
    i = directed_flag? nnodes*(nnodes-1) : (nnodes*(nnodes-1))/2;
    nw.duration_info.lasttoggle = (int *) calloc(i, sizeof(int));
  }
  else nw.duration_info.lasttoggle = NULL;

  /*Configure a Network*/
  nw.nnodes = nnodes;
  nw.nedges = 0; /* Edges will be added one by one */
  nw.directed_flag=directed_flag;
  nw.bipartite=bipartite;

  WtShuffleEdges(tails,heads,weights,nedges); /* shuffle to avoid worst-case performance */

  for(i = 0; i < nedges; i++) {
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

/*Takes vectors of doubles for edges; used only when constructing from inputparams. */

/* *** don't forget tail->head, so this function now accepts tails before heads */

WtNetwork WtNetworkInitializeD(double *tails, double *heads, double *weights, Edge nedges,
			     Vertex nnodes, int directed_flag, Vertex bipartite,
			     int lasttoggle_flag) {

  Vertex *itails=malloc(sizeof(Vertex)*nedges);
  Vertex *iheads=malloc(sizeof(Vertex)*nedges);
  
  for(Edge i=0; i<nedges; i++){
    itails[i]=tails[i];
    iheads[i]=heads[i];
  }

  WtNetwork nw=WtNetworkInitialize(itails,iheads,weights,nedges,nnodes,directed_flag,bipartite,lasttoggle_flag);

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

/*****************
 Edge EdgetreeSearch

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
 Edge EdgetreeMinimum

 Return the index of the WtTreeNode with the
 smallest value in the subtree rooted at x
*****************/
Edge WtEdgetreeMinimum (WtTreeNode *edges, Edge x) {
  Edge y;

  while ((y=(edges+x)->left) != 0)
    x=y;
  return x;
}


Edge DaveWtEdgetreeMinimum (WtTreeNode *edges, Edge x) {
  Edge y;
Rprintf("Here we go!\n");

Wtprintedge(x, edges);
  while ((y=(edges+x)->left) != 0) {
    x=y;
Wtprintedge(x, edges);
  }
  return x;
}

/*****************
 Edge ToggleEdge

 Toggle an edge:  Set it to the opposite of its current
 value.  Return 1 if edge added, 0 if deleted.
*****************/

/* *** don't forget tail->head, so this function now accepts tail before head */

int WtToggleEdge (Vertex tail, Vertex head, double weight, WtNetwork *nwp) 
{
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

  if (!(nwp->directed_flag) && tail > head) {
    Vertex temp;
    temp = tail; /*  Make sure tail<head always for undirected edges */
    tail = head;
    head = temp;
  }
  
  if(nwp->duration_info.lasttoggle){ /* Skip timestamps if no duration info. */
    if (nwp->directed_flag) 
      k = (head-1)*(nwp->nnodes-1) + tail - ((tail>head) ? 1:0) - 1; 
    else
      k = (head-1)*(head-2)/2 + tail - 1;    
    nwp->duration_info.lasttoggle[k] = nwp->duration_info.MCMCtimer;
  }

  if (WtAddEdgeToTrees(tail,head,weight,nwp))
    return 1;
  else 
    return 1 - WtDeleteEdgeFromTrees(tail,head,nwp);
}

/*****************
 long int ElapsedTime

 Return time since given (tail,head) was last toggled using
 ToggleEdgeWithTimestamp function
 *****************/

/* *** don't forget tail->head, so this function now accepts tail before head */

int WtElapsedTime (Vertex tail, Vertex head, WtNetwork *nwp) 
{
  Edge k;
  if (!(nwp->directed_flag) && tail > head) {
    Vertex temp;
    temp = tail; /*  Make sure tail<head always for undirected edges */
    tail = head;
    head = temp;
  }

  if(nwp->duration_info.lasttoggle){ /* Return 0 if no duration info. */
    if (nwp->directed_flag) 
      k = (head-1)*(nwp->nnodes-1) + tail - ((tail>head) ? 1:0) - 1; 
    else
      k = (head-1)*(head-2)/2 + tail - 1;    
    return nwp->duration_info.MCMCtimer - nwp->duration_info.lasttoggle[k];
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

void WtTouchEdge(Vertex tail, Vertex head, WtNetwork *nwp){
  unsigned int k;
  if(nwp->duration_info.lasttoggle){ /* Skip timestamps if no duration info. */
    if (nwp->directed_flag) 
      k = (head-1)*(nwp->nnodes-1) + tail - ((tail>head) ? 1:0) - 1; 
    else
      k = (head-1)*(head-2)/2 + tail - 1;    
    nwp->duration_info.lasttoggle[k] = nwp->duration_info.MCMCtimer;
  }
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

int WtAddEdgeToTrees(Vertex tail, Vertex head, double weight, WtNetwork *nwp){
  if (WtEdgetreeSearch(tail, head, nwp->outedges) == 0) {
    WtAddHalfedgeToTree(tail, head, weight, nwp->outedges, nwp->next_outedge);
    WtAddHalfedgeToTree(head, tail, weight, nwp->inedges, nwp->next_inedge);
    ++nwp->outdegree[tail];
    ++nwp->indegree[head];
    ++nwp->nedges;
    WtUpdateNextedge (nwp->inedges, &(nwp->next_inedge), nwp); 
    WtUpdateNextedge (nwp->outedges, &(nwp->next_outedge), nwp);
    return 1;
  }
  return 0;
}

/*****************
 void WtAddHalfedgeToTree:  Only called by WtAddEdgeToTrees
*****************/
void WtAddHalfedgeToTree (Vertex a, Vertex b, double weight, WtTreeNode *edges, Edge next_edge){
  WtTreeNode *eptr = edges+a, *newnode;
  Edge e;

  if (eptr->value==0) { /* This is the first edge for vertex a. */
    eptr->value=b;
    eptr->weight = weight;  /*  Add weight too */
    return;
  }
  (newnode = edges + next_edge)->value=b;  
  newnode->left = newnode->right = 0;
  newnode->weight=weight;  /*  Add weight too */
  /* Now find the parent of this new edge */
  for (e=a; e!=0; e=(b < (eptr=edges+e)->value) ? eptr->left : eptr->right);
  newnode->parent=eptr-edges;  /* Point from the new edge to the parent... */
  if (b < eptr->value)  /* ...and have the parent point back. */
    eptr->left=next_edge; 
  else
    eptr->right=next_edge;
}

/*****************
void UpdateNextedge
*****************/
void WtUpdateNextedge (WtTreeNode *edges, Edge *nextedge, WtNetwork *nwp) {
  int mult=2;
  
  while (++*nextedge < nwp->maxedges) {
    if (edges[*nextedge].value==0) return;
  }
  /* Reached end of allocated memory;  back to start and recheck for "holes" */
  for (*nextedge = (Edge)nwp->nnodes+1; *nextedge < nwp->maxedges; ++*nextedge) {
    if (edges[*nextedge].value==0) return;
  }
  /* There are no "holes" left, so this network overflows mem allocation */
  nwp->maxedges *= mult;
  nwp->inedges = (WtTreeNode *) realloc(nwp->inedges, 
					sizeof(WtTreeNode) * nwp->maxedges);
  memset(nwp->inedges+*nextedge,0,sizeof(WtTreeNode) * (nwp->maxedges-*nextedge));
  nwp->outedges = (WtTreeNode *) realloc(nwp->outedges, 
					 sizeof(WtTreeNode) * nwp->maxedges);
  memset(nwp->outedges+*nextedge,0,sizeof(WtTreeNode) * (nwp->maxedges-*nextedge));
}

/*****************
 int WtDeleteEdgeFromTrees

 Find and delete the edge from tail to head.  
 Return 1 if successful, 0 otherwise.  As with AddEdgeToTrees, this must
 be done once for outedges and once for inedges.
*****************/

/* *** don't forget tail->head, so this function now accepts tail before head */

int WtDeleteEdgeFromTrees(Vertex tail, Vertex head, WtNetwork *nwp){
  if (WtDeleteHalfedgeFromTree(tail, head, nwp->outedges,&(nwp->next_outedge))&&
      WtDeleteHalfedgeFromTree(head, tail, nwp->inedges, &(nwp->next_inedge))) {
    --nwp->outdegree[tail];
    --nwp->indegree[head];
    --nwp->nedges;
    return 1;
  }
  return 0;
}

/*****************
 int WtDeleteHalfedgeFromTree

 Delete the WtTreeNode with value b from the tree rooted at edges[a].
 Return 0 if no such WtTreeNode exists, 1 otherwise.  Also update the
 value of *next_edge appropriately.
*****************/
int WtDeleteHalfedgeFromTree(Vertex a, Vertex b, WtTreeNode *edges,
		     Edge *next_edge){ 
  Edge x, z, root=(Edge)a;
  WtTreeNode *xptr, *zptr, *ptr;

  if ((z=WtEdgetreeSearch(a, b, edges))==0)  /* z is the current WtTreeNode. */
    return 0; /* This edge doesn't exist, so return 0 */
  /* First, determine which node to splice out; this is z.  If the current
     z has two children, then we'll actually splice out its successor. */
  if ((zptr=edges+z)->left != 0 && zptr->right != 0) {
    z=WtEdgetreeSuccessor(edges, z);
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
  /* Clear z node, update *next_edge if necessary. */
  zptr->value=0;
  zptr->weight=0;
  if (z < *next_edge)
    *next_edge=z;
  return 1;
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
    Rprintf(" %d ",(edges+x)->value); 
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

/*****************
  int FindithEdge

  Find the ith edge in the Network *nwp and
  update the values of tail and head appropriately.  Return
  1 if successful, 0 otherwise.  
  Note that i is numbered from 1, not 0.  Thus, the maximum possible
  value of i is nwp->nedges.
******************/

/* *** don't forget tail->head, so this function now accepts tail before head */

int WtFindithEdge (Vertex *tail, Vertex *head, Edge i, WtNetwork *nwp) {
  Vertex taili=1;
  Edge e;

  if (i > nwp->nedges || i<=0)
    return 0;
  while (i > nwp->outdegree[taili]) {
    i -= nwp->outdegree[taili];
    taili++;
  }
  e=WtEdgetreeMinimum(nwp->outedges,taili);
  while (i-- > 1) {
    e=WtEdgetreeSuccessor(nwp->outedges, e);
  }
  *tail = taili;
  *head = nwp->outedges[e].value;
  return 1;
}

/*****************
  int WtFindithnonEdge

  Find the ith nonedge in the WtNetwork *nwp and
  update the values of tail and head appropriately.  Return
  1 if successful, 0 otherwise.  
  Note that i is numbered from 1, not 0.  Thus, the maximum possible
  value of i is (ndyads - nwp->nedges).
******************/

/* This function is not yet written.  It's not clear whether it'll
   be needed. */      
  /* *** but if it is needed, don't forget,  tail -> head */

/* int WtFindithnonEdge (Vertex *tail, Vertex *head, Edge i, WtNetwork *nwp) {
} */


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
	if(weights) weights[nextedge] = EdgeWeight(tails[nextedge],heads[nextedge],nwp);
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
	  if(weights) weights[nextedge] = EdgeWeight(k,v,nwp);
          nextedge++;
        }else{
          tails[nextedge] = v;
          heads[nextedge] = k;
	  if(weights) weights[nextedge] = EdgeWeight(v,k,nwp);
          nextedge++;
        }
      }
    }
  }
  return nextedge;
}

void WtShuffleEdges(Vertex *tails, Vertex *heads, double *weights, Edge nedges){
  /* *** don't forget,  tail -> head */
  for(Edge i = nedges; i > 0; i--) {
    Edge j = i * unif_rand();  /* shuffle to avoid worst-case performance */
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

/*****************
 long int EdgeWeight

 Return weight of (tail,head) in a WtNetwork
 *****************/
double EdgeWeight (Vertex tail, Vertex head, WtNetwork *nwp) 
{
  Edge k;
  if (!(nwp->directed_flag) && tail > head) {
    Vertex temp;
    temp = tail; /*  Make sure tail<head always for undirected edges */
    tail = head;
    head = temp;
  }

  if (nwp->directed_flag) 
    k = (head-1)*(nwp->nnodes-1) + tail - ((tail>head) ? 1:0) - 1; 
  else
    k = (head-1)*(head-2)/2 + tail - 1;    
  return nwp->outedges[k].weight;
}
