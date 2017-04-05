#include "edgetree.h"

/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */


/*******************
 Network NetworkInitialize

 Initialize, construct binary tree version of network.  Note
 that the 0th TreeNode in the array is unused and should 
 have all its values set to zero
*******************/
Network NetworkInitialize(Vertex *tails, Vertex *heads, Edge nedges, 
			  Vertex nnodes, int directed_flag, Vertex bipartite,
			  int lasttoggle_flag) {

  /* *** don't forget, tail -> head */

  Network nw;

  nw.next_inedge = nw.next_outedge = (Edge)nnodes+1;
  /* Calloc will zero the allocated memory for us, probably a lot
     faster. */
  nw.outdegree = (Vertex *) calloc((nnodes+1),sizeof(Vertex));
  nw.indegree  = (Vertex *) calloc((nnodes+1),sizeof(Vertex));
  nw.maxedges = MAX(nedges,1)+nnodes+2; /* Maybe larger than needed? */
  nw.inedges = (TreeNode *) calloc(nw.maxedges,sizeof(TreeNode));
  nw.outedges = (TreeNode *) calloc(nw.maxedges,sizeof(TreeNode));

  GetRNGstate();  /* R function enabling uniform RNG */

  if(lasttoggle_flag){
    nw.duration_info.MCMCtimer=0;
    nw.duration_info.lasttoggle = (int *) calloc(directed_flag? nnodes*(nnodes-1) : (nnodes*(nnodes-1))/2, sizeof(int));
  }
  else nw.duration_info.lasttoggle = NULL;

  /*Configure a Network*/
  nw.nnodes = nnodes;
  nw.nedges = 0; /* Edges will be added one by one */
  nw.directed_flag=directed_flag;
  nw.bipartite=bipartite;

  ShuffleEdges(tails,heads,nedges); /* shuffle to avoid worst-case performance */

  for(Edge i = 0; i < nedges; i++) {
    Vertex tail=tails[i], head=heads[i];
    if (!directed_flag && tail > head) 
      AddEdgeToTrees(head,tail,&nw); /* Undir edges always have tail < head */ 
    else 
      AddEdgeToTrees(tail,head,&nw);
  }
  PutRNGstate();
  return nw;
}


/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */

/*Takes vectors of doubles for edges; used only when constructing from inputparams. */
Network NetworkInitializeD(double *tails, double *heads, Edge nedges,
			  Vertex nnodes, int directed_flag, Vertex bipartite,
			  int lasttoggle_flag) {

  /* *** don't forget, tail -> head */

  Vertex *itails=(Vertex*)malloc(sizeof(Vertex)*nedges);
  Vertex *iheads=(Vertex*)malloc(sizeof(Vertex)*nedges);
  
  for(Edge i=0; i<nedges; i++){
    itails[i]=tails[i];
    iheads[i]=heads[i];
  }

  Network nw=NetworkInitialize(itails,iheads,nedges,nnodes,directed_flag,bipartite,lasttoggle_flag);

  free(itails);
  free(iheads);
  return nw;
}

/*******************
 void NetworkDestroy
*******************/
void NetworkDestroy(Network *nwp) {
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
 Network NetworkCopy
*****************/
Network *NetworkCopy(Network *dest, Network *src){
  Vertex nnodes = dest->nnodes = src->nnodes;
  dest->next_inedge = src->next_inedge;
  dest->next_outedge = src->next_outedge;

  dest->outdegree = (Vertex *) malloc((nnodes+1)*sizeof(Vertex));
  memcpy(dest->outdegree, src->outdegree, (nnodes+1)*sizeof(Vertex));
  dest->indegree = (Vertex *) malloc((nnodes+1)*sizeof(Vertex));
  memcpy(dest->indegree, src->indegree, (nnodes+1)*sizeof(Vertex));

  Vertex maxedges = dest->maxedges = src->maxedges;

  dest->inedges = (TreeNode *) malloc(maxedges*sizeof(TreeNode));
  memcpy(dest->inedges, src->inedges, maxedges*sizeof(TreeNode));
  dest->outedges = (TreeNode *) malloc(maxedges*sizeof(TreeNode));
  memcpy(dest->outedges, src->outedges, maxedges*sizeof(TreeNode));

  int directed_flag = dest->directed_flag = src->directed_flag;

  if(src->duration_info.lasttoggle){
    dest->duration_info.MCMCtimer=src->duration_info.MCMCtimer;
    dest->duration_info.lasttoggle = (int *) calloc(directed_flag? nnodes*(nnodes-1) : (nnodes*(nnodes-1))/2, sizeof(int));
    memcpy(dest->duration_info.lasttoggle, src->duration_info.lasttoggle,(directed_flag? nnodes*(nnodes-1) : (nnodes*(nnodes-1))/2) * sizeof(int));
  }
  else dest->duration_info.lasttoggle = NULL;

  dest->nedges = src->nedges;
  dest->bipartite = src->bipartite;

  return dest;
}

/*****************
 Edge EdgetreeSearch

 Check to see if there's a TreeNode with value b 
 in the tree rooted at edges[a].  Return i such that 
 edges[i] is that TreeNode, or 0 if none.
*****************/
Edge EdgetreeSearch (Vertex a, Vertex b, TreeNode *edges) {
  TreeNode *es;
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
Edge EdgetreeSuccessor (TreeNode *edges, Edge x) {
  TreeNode *ptr;
  Edge y;

  if ((y=(ptr=edges+x)->right) != 0) 
    return EdgetreeMinimum (edges, y);
  while ((y=ptr->parent)!=0 && x==(ptr=edges+y)->right) 
    x=y;
  return y; 
}   

/*****************
 Edge EdgetreeMinimum

 Return the index of the TreeNode with the
 smallest value in the subtree rooted at x
*****************/
Edge EdgetreeMinimum (TreeNode *edges, Edge x) {
  Edge y;

  while ((y=(edges+x)->left) != 0)
    x=y;
  return x;
}



/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */



/*****************
 Edge ToggleEdge

 Toggle an edge:  Set it to the opposite of its current
 value.  Return 1 if edge added, 0 if deleted.
*****************/
int ToggleEdge (Vertex tail, Vertex head, Network *nwp) 
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
int ToggleEdgeWithTimestamp(Vertex tail, Vertex head, Network *nwp){
  Edge k;

  /* don't forget, tails < heads in undirected networks now  */
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
  
  if (AddEdgeToTrees(tail,head,nwp))
    return 1;
  else 
    return 1 - DeleteEdgeFromTrees(tail,head,nwp);
}



/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */



/*****************
 long int ElapsedTime

 Return time since given (tail,head) was last toggled using
 ToggleEdgeWithTimestamp function
*****************/
int ElapsedTime(Vertex tail, Vertex head, Network *nwp){
  Edge k;
  /* don't forget, tails < heads now in undirected networks */
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



/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */


/*****************
 void TouchEdge

 Named after the UNIX "touch" command.
 Set an edge's time-stamp to the current MCMC time.
 *****************/

void TouchEdge(Vertex tail, Vertex head, Network *nwp){
  unsigned int k;
  if(nwp->duration_info.lasttoggle){ /* Skip timestamps if no duration info. */
    if (nwp->directed_flag) 
      k = (head-1)*(nwp->nnodes-1) + tail - ((tail>head) ? 1:0) - 1; 
    else
      k = (head-1)*(head-2)/2 + tail - 1;    
    nwp->duration_info.lasttoggle[k] = nwp->duration_info.MCMCtimer;
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
int AddEdgeToTrees(Vertex tail, Vertex head, Network *nwp){
  if (EdgetreeSearch(tail, head, nwp->outedges) == 0) {
    AddHalfedgeToTree(tail, head, nwp->outedges, nwp->next_outedge);
    AddHalfedgeToTree(head, tail, nwp->inedges, nwp->next_inedge);
    ++nwp->outdegree[tail];
    ++nwp->indegree[head];
    ++nwp->nedges;
    UpdateNextedge (nwp->inedges, &(nwp->next_inedge), nwp); 
    UpdateNextedge (nwp->outedges, &(nwp->next_outedge), nwp);
    return 1;
  }
  return 0;
}

/*****************
 void AddHalfedgeToTree:  Only called by AddEdgeToTrees
*****************/
void AddHalfedgeToTree (Vertex a, Vertex b, TreeNode *edges, Edge next_edge){
  TreeNode *eptr = edges+a, *newnode;
  Edge e;

  if (eptr->value==0) { /* This is the first edge for vertex a. */
    eptr->value=b;
    return;
  }
  (newnode = edges + next_edge)->value=b;  
  newnode->left = newnode->right = 0;
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
void UpdateNextedge (TreeNode *edges, Edge *nextedge, Network *nwp) {
  int mult=2;
  /*TreeNode *tmp_in, *tmp_out; */
  
  while (++*nextedge < nwp->maxedges) {
    if (edges[*nextedge].value==0) return;
  }
  /* Reached end of allocated memory;  back to start and recheck for "holes" */
  for (*nextedge = (Edge)nwp->nnodes+1; *nextedge < nwp->maxedges; ++*nextedge) {
    if (edges[*nextedge].value==0) return;
  }
  /* There are no "holes" left, so this network overflows mem allocation */
  nwp->maxedges *= mult;
  nwp->inedges = (TreeNode *) realloc(nwp->inedges, 
                                      sizeof(TreeNode) * nwp->maxedges);
  memset(nwp->inedges+*nextedge,0,sizeof(TreeNode) * (nwp->maxedges-*nextedge));
  nwp->outedges = (TreeNode *) realloc(nwp->outedges, 
                                       sizeof(TreeNode) * nwp->maxedges);
  memset(nwp->outedges+*nextedge,0,sizeof(TreeNode) * (nwp->maxedges-*nextedge));
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
int DeleteEdgeFromTrees(Vertex tail, Vertex head, Network *nwp){
  if (DeleteHalfedgeFromTree(tail, head, nwp->outedges,&(nwp->next_outedge))&&
      DeleteHalfedgeFromTree(head, tail, nwp->inedges, &(nwp->next_inedge))) {
    --nwp->outdegree[tail];
    --nwp->indegree[head];
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
int DeleteHalfedgeFromTree(Vertex a, Vertex b, TreeNode *edges,
		     Edge *next_edge){
  Edge x, z, root=(Edge)a;
  TreeNode *xptr, *zptr, *ptr;

  if ((z=EdgetreeSearch(a, b, edges))==0)  /* z is the current TreeNode. */
    return 0;  /* This edge doesn't exist, so return 0 */
  /* First, determine which node to splice out; this is z.  If the current
     z has two children, then we'll actually splice out its successor. */
  if ((zptr=edges+z)->left != 0 && zptr->right != 0) {
    z=EdgetreeSuccessor(edges, z);  
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
  /* Clear z node, update *next_edge if necessary. */
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
  value of i is nwp->nedges.
******************/
int FindithEdge (Vertex *tail, Vertex *head, Edge i, Network *nwp) {
  Vertex head1=1;
  Edge e;

  if (i > nwp->nedges || i<=0)
    return 0;
  while (i > nwp->outdegree[head1]) {
    i -= nwp->outdegree[head1];
    head1++;
  }
  e=EdgetreeMinimum(nwp->outedges,head1);
  while (i-- > 1) {
    e=EdgetreeSuccessor(nwp->outedges, e);
  }
  *tail = head1;
  *head = nwp->outedges[e].value;
  return 1;
}

/*****************
  int FindithnonEdge

  Find the ith nonedge in the Network *nwp and
  update the values of tail and head appropriately.  Return
  1 if successful, 0 otherwise.  
  Note that i is numbered from 1, not 0.  Thus, the maximum possible
  value of i is (ndyads - nwp->nedges).
******************/

/* This function is not yet written.  It's not clear whether it'll
   be needed. */      

/* int FindithnonEdge (Vertex *tail, Vertex *head, Edge i, Network *nwp) {
} */




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
  if (nwp->directed_flag) {
    for (Vertex v=1; v<=nwp->nnodes; v++){
      for(Vertex e = EdgetreeMinimum(nwp->outedges, v);
      nwp->outedges[e].value != 0 && nextedge < nmax;
      e = EdgetreeSuccessor(nwp->outedges, e)){
        tails[nextedge] = v;
        heads[nextedge] = nwp->outedges[e].value;
        nextedge++;
      }
    }
  }else{
    for (Vertex v=1; v<=nwp->nnodes; v++){
      for(Vertex e = EdgetreeMinimum(nwp->outedges, v);
      nwp->outedges[e].value != 0 && nextedge < nmax;
      e = EdgetreeSuccessor(nwp->outedges, e)){
        Vertex k = nwp->outedges[e].value;
        if(v < k){
          tails[nextedge] = k;
          heads[nextedge] = v;
          nextedge++;
        }else{
          tails[nextedge] = v;
          heads[nextedge] = k;
          nextedge++;
        }
      }
    }
  }
  return nextedge;
}


/* *** don't forget, edges are now given by tails -> heads, and as
       such, the function definitions now require tails to be passed
       in before heads */


void ShuffleEdges(Vertex *tails, Vertex *heads, Edge nedges){
  for(Edge i = nedges; i > 0; i--) {
    Edge j = (double) i * unif_rand();  /* shuffle to avoid worst-case performance */
    Vertex tail = tails[j];
    Vertex head = heads[j];
    tails[j] = tails[i-1];
    heads[j] = heads[i-1];
    tails[i-1] = tail;
    heads[i-1] = head;
  }
}
