/*  File inst/include/ergm_edgetree.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2019 Statnet Commons
 */
#ifndef _ERGM_EDGETREE_H_
#define _ERGM_EDGETREE_H_

#include "ergm_edgetree_common.do_not_include_directly.h"

/* Ensure that tail < head for undriected networks. */
#define ENSURE_TH_ORDER							\
  if(!(nwp->directed_flag) && tail>head){				\
    Vertex temp;							\
    temp = tail;							\
    tail = head;							\
    head = temp;							\
  }

/*  TreeNode is a binary tree structure, which is how the edgelists 
    are stored.  The root of the tree for vertex i will be inedges[i]
    or outedges[i].  inedges[0] and outedges[0] are unused, since the
    index 0 will indicate no link.  Indices are long unsigned integers,
    which means networks can contain 2^32-1= 4294967295 edges (enough to
    hold all edges in a 92682-vertex undirected network or a 65536-vertex
    directed network, assuming no multiple edges or loops), though for this
    MAXEDGES must be adjusted accordingly.
*/
typedef struct TreeNodestruct {
  Vertex value;      /*  the vertex at the other end of this edge  */
  Edge parent;   /*  parent of this node in the tree (0 for root) */
  Edge left;     /*  left child (0 if none)  */
  Edge right;    /*  right child (0 if none) */
} TreeNode;

/* Network is a structure containing all essential elements
   of a given network; it is a slightly rewritten version of the old Gptr,
   with some changes of awkard things, deletion of unnecessary things, and
   a new name more reflective of what it does!

   Some of the fields in a Network structure are:
   inedges and outedges are arrays of TreeNode that are used to 
     store all of the incoming and outgoing edges, respectively. 
   directed_flag is 1 or 0, depending on whether or not the 
     network is directed. 
   last_inedge and last_outedge are continually updated to give
     the highest index of an edge object being used.  
   outdegree[] and indegree[] are continually updated to give
     the appropriate degree values for each vertex.  These should
     point to Vertex-vectors of length nnodes+1.  
   value:  optional value(s) associated with this network 
*/
typedef struct Networkstruct {
  TreeNode *inedges;
  TreeNode *outedges;
  int directed_flag;
  Vertex bipartite;  
  Vertex nnodes;
  Edge nedges;
  Edge last_inedge;
  Edge last_outedge;
  Vertex *indegree;
  Vertex *outdegree;
  double *value;  
  Dur_Inf duration_info;
  Edge maxedges;
} Network;


/* *** don't forget,  tails -> heads, so all the functions below using
   heads & tails, now list tails before heads */

/* Initialization and destruction. */
Network *NetworkInitialize(Vertex *tails, Vertex *heads, Edge nedges,
			   Vertex nnodes, int directed_flag, Vertex bipartite,
			   int lasttoggle_flag, int time, int *lasttoggle);
void NetworkDestroy(Network *nwp);
Network *NetworkInitializeD(double *tails, double *heads, Edge nedges,
			    Vertex nnodes, int directed_flag, Vertex bipartite,
			    int lasttoggle_flag, int time, int *lasttoggle);

Network *NetworkCopy(Network *src);

/* /\* Accessors. *\/ */
/* static inline Edge EdgetreeSearch (Vertex a, Vertex b, TreeNode *edges); */
/* static inline Edge EdgetreeSuccessor (TreeNode *edges, Edge x); */
/* static inline Edge EdgetreePredecessor (TreeNode *edges, Edge x); */
/* static inline Edge EdgetreeMinimum (TreeNode *edges, Edge x); */
/* static inline Edge EdgetreeMaximum (TreeNode *edges, Edge x); */

/* Modifiers. */

/* *** don't forget,  tails -> heads, so all the functions below using
   heads & tails, now list tails before heads */

void SetEdge (Vertex tail, Vertex head, unsigned int weight, Network *nwp);
void SetEdgeWithTimestamp (Vertex tail, Vertex head, unsigned int weight, Network *nwp);
int ToggleEdge (Vertex tail, Vertex head, Network *nwp);
int ToggleEdgeWithTimestamp (Vertex tail, Vertex head, Network *nwp);
int AddEdgeToTrees(Vertex tail, Vertex head, Network *nwp);
void AddHalfedgeToTree (Vertex a, Vertex b, TreeNode *edges, Edge *last_edge);
void CheckEdgetreeFull (Network *nwp);
int DeleteEdgeFromTrees(Vertex tail, Vertex head, Network *nwp);
int DeleteHalfedgeFromTree(Vertex a, Vertex b, TreeNode *edges,
		     Edge *last_edge);
void RelocateHalfedge(Edge from, Edge to, TreeNode *edges);

/* /\* Duration functions. *\/ */
/* static inline int ElapsedTime(Vertex tail, Vertex head, Network *nwp); */
void TouchEdge(Vertex tail, Vertex head, Network *nwp);

#include "ergm_edgetree_inline.do_not_include_directly.h"

/* Utility functions. */
int FindithEdge(Vertex *tail, Vertex *head, Edge i, Network *nwp);
int GetRandEdge(Vertex *tail, Vertex *head, Network *nwp);
int FindithNondge(Vertex *tail, Vertex *head, Dyad i, Network *nwp);
int GetRandNonedge(Vertex *tail, Vertex *head, Network *nwp);
void printedge(Edge e, TreeNode *edges);
void InOrderTreeWalk(TreeNode *edges, Edge x);
void NetworkEdgeList(Network *nwp);
void ShuffleEdges(Vertex *tails, Vertex *heads, Edge nedges);
void DetShuffleEdges(Vertex *tails, Vertex *heads, Edge nedges);

/* Others... */
Edge DesignMissing (Vertex a, Vertex b, Network *mnwp);
Edge EdgeTree2EdgeList(Vertex *tails, Vertex *heads, Network *nwp, Edge nmax);

#endif
