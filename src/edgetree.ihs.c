#include "edgeTree.ihs.h"

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
 Network NetworkInitialize

 Initialize, construct binary tree version of network.  Note
 that the 0th TreeNode in the array is unused and should 
 have all its values set to zero
*******************/
WtNetwork WtNetworkInitialize(int *heads, int *tails, double *weights,
			      int nedges, int nnodes, int directed_flag, int bipartite) {
  Edge i, j;
  Vertex h, t;
  double w;
  WtNetwork nw;

  nw.next_inedge = nw.next_outedge = (Edge)nnodes+1;
  nw.outdegree = (Vertex *) malloc(sizeof(Vertex) * (nnodes+1));
  nw.indegree  = (Vertex *) malloc(sizeof(Vertex) * (nnodes+1));
  nw.inedges = (WtTreeNode *) malloc(sizeof(WtTreeNode) * MAXEDGES);
  nw.outedges = (WtTreeNode *) malloc(sizeof(WtTreeNode) * MAXEDGES);
  
  nw.duration_info.MCMCtimer=0;
  i = directed_flag? nnodes*(nnodes-1) : (nnodes*(nnodes-1))/2;
  nw.duration_info.lasttoggle = (int *) malloc(sizeof(int) * i);

  for (i=0; i<=nnodes; i++) {
    nw.inedges[i].value = nw.outedges[i].value = 0;
    nw.inedges[i].parent = nw.outedges[i].parent = 0;
    nw.inedges[i].left = nw.outedges[i].left = 0;
    nw.inedges[i].right = nw.outedges[i].right = 0;
    nw.inedges[i].weight = nw.outedges[i].weight = 0.0;
    nw.indegree[i] = nw.outdegree[i] = 0;
  }
  
  for (; i<MAXEDGES; i++)
    nw.inedges[i].value = nw.outedges[i].value = 0;

  /*Configure a Network*/
  nw.nnodes = nnodes;
  nw.nedges = 0; /* Edges will be added one by one */
  nw.directed_flag=directed_flag;
  nw.bipartite=bipartite;

  for(i = nedges; i > 0; i--) {
    j = i * unif_rand();  /* shuffle to avoid worst-case performance */
    h = (Vertex)heads[j];
    t = (Vertex)tails[j];
    w = weights[j];
    heads[j] = heads[i-1];
    tails[j] = tails[i-1];
    weights[j] = weights[i-1];
    heads[i-1] = h;
    tails[i-1] = t;
    weights[i-1] = w;
    if (!directed_flag && h > t) 
      WtAddEdgeToTrees(t,h,w,&nw); /* Undir edges always have head < tail */ 
    else 
      WtAddEdgeToTrees(h,t,w,&nw);
  }
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
  free (nwp->duration_info.lasttoggle);
}

/*****************
 Edge EdgetreeSearch

 Check to see if there's a TreeNode with value b 
 in the tree rooted at edges[a].  Return i such that 
 edges[i] is that TreeNode, or 0 if none.
*****************/
Edge WtEdgetreeSearch (Vertex a, Vertex b, WtTreeNode *edges) {
  WtTreeNode *es;
  Edge e = a;
  Vertex v;

  es = edges + e;
  v = es->value;
  while(e != 0 && b != v)
    {
      e = (b<v)?  es->left : es->right;
      es = edges + e;
      v = es->value;
    }
  return e;
}

/*****************
 Edge EdgetreeSuccessor

 Return the index of the TreeNode with the smallest value
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

 Return the index of the TreeNode with the
 smallest value in the subtree rooted at x
*****************/
Edge WtEdgetreeMinimum (WtTreeNode *edges, Edge x) {
  Edge y;

  while ((y=(edges+x)->left) != 0)
    x=y;
  return x;
}

/*****************
 Edge ToggleEdge

 Toggle an edge:  Set it to the opposite of its current
 value.  Return 1 if edge added, 0 if deleted.
*****************/
int WtToggleEdge (Vertex head, Vertex tail, double weight, WtNetwork *nwp) 
{
  if (!(nwp->directed_flag) && head > tail) {
    Vertex temp;
    temp = head; /*  Make sure head<tail always for undirected edges */
    head = tail;
    tail = temp;
  }
  if (WtAddEdgeToTrees(head,tail,weight,nwp))
    return 1;
  else 
    return 1 - WtDeleteEdgeFromTrees(head,tail,nwp);
}


/*****************
 Edge ToggleEdgeWithTimestamp
 By MSH 11/26/06

 Same as ToggleEdge, but this time with the additional
 step of updating the matrix of 'lasttoggle' times
 *****************/
int WtToggleEdgeWithTimestamp (Vertex head, Vertex tail, double weight, WtNetwork *nwp) 
{
  Edge k;
  if (!(nwp->directed_flag) && head > tail) {
    Vertex temp;
    temp = head; /*  Make sure head<tail always for undirected edges */
    head = tail;
    tail = temp;
  }
  if(nwp->duration_info.lasttoggle){ /* Skip timestamps if no duration info. */
    if (nwp->directed_flag) 
      k = (tail-1)*(nwp->nnodes-1) + head - ((head>tail) ? 1:0) - 1; 
    else
      k = (tail-1)*(tail-2)/2 + head - 1;    
    nwp->duration_info.lasttoggle[k] = nwp->duration_info.MCMCtimer;
  }

  if (WtAddEdgeToTrees(head,tail,weight,nwp))
    return 1;
  else 
    return 1 - WtDeleteEdgeFromTrees(head,tail,nwp);
}

/*****************
 long int ElapsedTime

 Return time since given (head,tail) was last toggled using
 ToggleEdgeWithTimestamp function
 *****************/
int ElapsedTime (Vertex head, Vertex tail, Network *nwp) 
{
  Edge k;
  if (!(nwp->directed_flag) && head > tail) {
    Vertex temp;
    temp = head; /*  Make sure head<tail always for undirected edges */
    head = tail;
    tail = temp;
  }

  if(nwp->duration_info.lasttoggle){ /* Return 0 if no duration info. */
    if (nwp->directed_flag) 
      k = (tail-1)*(nwp->nnodes-1) + head - ((head>tail) ? 1:0) - 1; 
    else
      k = (tail-1)*(tail-2)/2 + head - 1;    
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

void TouchEdge(Vertex head, Vertex tail, Network *nwp){
  unsigned int k;
  if(nwp->duration_info.lasttoggle){ /* Skip timestamps if no duration info. */
    if (nwp->directed_flag) 
      k = (tail-1)*(nwp->nnodes-1) + head - ((head>tail) ? 1:0) - 1; 
    else
      k = (tail-1)*(tail-2)/2 + head - 1;    
    nwp->duration_info.lasttoggle[k] = nwp->duration_info.MCMCtimer;
  }
}

/*****************
 Edge AddEdgeToTrees

 Add an edge from head to tail after checking to see
 if it's legal. Return 1 if edge added, 0 otherwise.  Since each
 "edge" should be added to both the list of outedges and the list of 
 inedges, this actually involves two calls to AddHalfEdgeToTree (hence
 "Trees" instead of "Tree" in the name of this function).
*****************/
int WtAddEdgeToTrees(Vertex head, Vertex tail, double weight, WtNetwork *nwp){
  if (WtEdgetreeSearch(head, tail, nwp->outedges) == 0) {
    WtAddHalfedgeToTree(head, tail, weight, nwp->outedges, 
			&(nwp->next_outedge));
    WtAddHalfedgeToTree(tail, head, weight, nwp->inedges, &(nwp->next_inedge));
    ++nwp->outdegree[head];
    ++nwp->indegree[tail];
    ++nwp->nedges;
    return 1;
  }
  return 0;
}

/*****************
 Edge AddHalfedgeToTree
*****************/
void WtAddHalfedgeToTree (Vertex a, Vertex b, double weight, 
			  WtTreeNode *edges, Edge *next_edge) 
{  /*  See comments in AddHalfedgeToTree.  */
  WtTreeNode *eptr = edges+a, *newnode;
  Edge e;
  
  if (eptr->value==0) { 
    eptr->value=b;
    eptr->weight = weight;  /*  Add weight too */
    return;
  }
  (newnode=edges+*next_edge)->value=b;  
  newnode->left=newnode->right=0;
  newnode->weight=weight;  /*  Add weight too */
  for (e=a; e!=0; e=(b < (eptr=edges+e)->value) ? eptr->left : eptr->right);
  newnode->parent=eptr-edges;
  if (b < eptr->value)
    eptr->left=*next_edge;
  else
    eptr->right=*next_edge;
  while (++*next_edge<MAXEDGES && edges[*next_edge].value!=0);
}

/*****************
 int DeleteEdgeFromTrees

 Find and delete the edge from head to tail.  
 Return 1 if successful, 0 otherwise.  As with AddEdgeToTrees, this must
 be done once for outedges and once for inedges.
*****************/
int WtDeleteEdgeFromTrees(Vertex head, Vertex tail, WtNetwork *nwp){
  if (WtDeleteHalfedgeFromTree(head, tail, nwp->outedges,
			       &(nwp->next_outedge))&&
      WtDeleteHalfedgeFromTree(tail, head, nwp->inedges, 
			       &(nwp->next_inedge))) {
    --nwp->outdegree[head];
    --nwp->indegree[tail];
    --nwp->nedges;
    return 1;
  }
  return 0;
}

/*****************
 int DeleteHalfedgeFromTree

 Delete the TreeNode with value b from the tree rooted at edges[a].
 Return 0 if no such TreeNode exists, 1 otherwise.  Also update the
 value of *next_edge appropriately.
*****************/
int WtDeleteHalfedgeFromTree(Vertex a, Vertex b, WtTreeNode *edges,
		     Edge *next_edge){ 
  /* see comments in DeleteHalfedgeFromTree */
  Edge x, z, root=(Edge)a;
  WtTreeNode *xptr, *zptr, *ptr;

  if ((z=WtEdgetreeSearch(a, b, edges))==0)  /* z is the current TreeNode. */
    return 0;
  if ((zptr=edges+z)->left != 0 && zptr->right != 0) {
    z=WtEdgetreeSuccessor(edges, z);
    zptr->value = (ptr=edges+z)->value;
    zptr=ptr;
  }
  if ((x=zptr->left) == 0)
    x = zptr->right;
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
  zptr->value=0;
  if (z < *next_edge)
    *next_edge=z;
  return 1;
}

/*****************
 void printedge

 Diagnostic routine that prints out the contents
 of the specified TreeNode (used for debugging).  
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
 void NetworkEdgeList

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
 void InOrderTreeWalk

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
  update the values of head and tail appropriately.  Return
  1 if successful, 0 otherwise.  
  Note that i is numbered from 1, not 0.  Thus, the maximum possible
  value of i is nwp->nedges.
******************/
int WtFindithEdge (Vertex *head, Vertex *tail, Edge i, WtNetwork *nwp) {
  Vertex h=1;
  Edge e;

  if (i > nwp->nedges || i<=0)
    return 0;
  while (i > nwp->outdegree[h]) {
    i -= nwp->outdegree[h];
    h++;
  }
  e=WtEdgetreeMinimum(nwp->outedges,h);
  while (i-- > 1) {
    e=WtEdgetreeSuccessor(nwp->outedges, e);
  }
  *head = h;                             
  *tail = nwp->outedges[e].value;
  return 1;
}


/*****************
 long int EdgeWeight

 Return weight of (head,tail) in a WtNetwork
 *****************/
double EdgeWeight (Vertex head, Vertex tail, WtNetwork *nwp) 
{
  Edge k;
  if (!(nwp->directed_flag) && head > tail) {
    Vertex temp;
    temp = head; /*  Make sure head<tail always for undirected edges */
    head = tail;
    tail = temp;
  }

  if (nwp->directed_flag) 
    k = (tail-1)*(nwp->nnodes-1) + head - ((head>tail) ? 1:0) - 1; 
  else
    k = (tail-1)*(tail-2)/2 + head - 1;    
  return nwp->outedges[k].weight;
}


