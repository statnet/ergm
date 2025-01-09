/*  File inst/include/ergm_wtedgetree.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#ifndef _ERGM_WTEDGETREE_H_
#define _ERGM_WTEDGETREE_H_

#include "ergm_edgetree_common.do_not_include_directly.h"

/* WtTreeNode is just like TreeNode but with an extra field for a
   weight, or value, that might be associated with the node */
typedef struct WtTreeNodestruct {
  Vertex value;      /*  the vertex at the other end of the edge  */
  Edge parent;   /*  parent of this node in the tree (0 for root) */
  Edge left;     /*  left child (0 if none)  */
  Edge right;    /*  right child (0 if none) */
  double weight;
} WtTreeNode;

/* WtNetwork is a structure just like Network except it is for a network with 
   weighted (valued) edges.  */
typedef struct WtNetworkstruct {
  WtTreeNode *inedges;
  WtTreeNode *outedges;
  int directed_flag;
  Vertex bipartite;  
  Vertex nnodes;
  Edge nedges;
  Edge last_inedge;
  Edge last_outedge;
  Vertex *indegree;
  Vertex *outdegree;
  const char *eattrname;
  Edge maxedges;

  unsigned int n_on_edge_change;
  unsigned int max_on_edge_change;
  void (**on_edge_change)(Vertex, Vertex, double, void*, struct WtNetworkstruct*, double);
  void **on_edge_change_payload;
} WtNetwork;
typedef void (*OnWtNetworkEdgeChange)(Vertex, Vertex, double, void*, WtNetwork*, double);

/* Initialization and destruction. */
WtNetwork *WtNetworkInitialize(Vertex *tails, Vertex *heads, double *weights, Edge nedges,
			       Vertex nnodes, int directed_flag, Vertex bipartite,
			       int lasttoggle_flag, int time, int *lasttoggle);
void WtNetworkDestroy(WtNetwork *nwp);
WtNetwork *WtNetworkInitializeD(double *tails, double *heads, double *weights, Edge nedges,
				Vertex nnodes, int directed_flag, Vertex bipartite,
				int lasttoggle_flag, int time, int *lasttoggle);

WtNetwork *WtNetworkCopy(WtNetwork *src);

SEXP WtNetwork2Redgelist(WtNetwork *nwp);
WtNetwork *Redgelist2WtNetwork(SEXP elR, Rboolean empty);

/* /\* Accessors. *\/ */
/* static inline Edge WtEdgetreeSearch (Vertex a, Vertex b, WtTreeNode *edges); */
/* static inline double WtGetEdge (Vertex tail, Vertex head, WtNetwork *nwp); */
/* static inline Edge WtEdgetreeSuccessor (WtTreeNode *edges, Edge x); */
/* static inline Edge WtEdgetreePredecessor (WtTreeNode *edges, Edge x); */
/* static inline Edge WtEdgetreeMinimum (WtTreeNode *edges, Edge x); */
/* static inline Edge WtEdgetreeMaximum (WtTreeNode *edges, Edge x); */

/* Modifiers. */

/* *** don't forget,  tails -> heads, so all the functions below using
   heads & tails, now list tails before heads */

void WtSetEdge (Vertex tail, Vertex head, double weight, WtNetwork *nwp);
int WtToggleEdge (Vertex tail, Vertex head, double weight, WtNetwork *nwp);
void WtAddEdgeToTrees(Vertex tail, Vertex head, double weight, WtNetwork *nwp);
/* void WtAddHalfedgeToTree (Vertex a, Vertex b, double weight, WtTreeNode *edges, Edge *last_edge); */
/* void WtCheckEdgetreeFull (WtNetwork *nwp); */
int WtDeleteEdgeFromTrees(Vertex tail, Vertex head, WtNetwork *nwp);
/* int WtDeleteHalfedgeFromTree(Vertex a, Vertex b, WtTreeNode *edges, */
/* 		     Edge *last_edge); */
/* void WtRelocateHalfedge(Edge from, Edge to, WtTreeNode *edges); */

/* Callback management. */
void AddOnWtNetworkEdgeChange(WtNetwork *nwp, OnWtNetworkEdgeChange callback, void *payload, unsigned int pos);
void DeleteOnWtNetworkEdgeChange(WtNetwork *nwp, OnWtNetworkEdgeChange callback, void *payload);

#include "ergm_wtedgetree_inline.do_not_include_directly.h"

/* Utility functions. */
int WtFindithEdge (Vertex *tail, Vertex *head, double *weight, Edge i, WtNetwork *nwp);
int WtGetRandEdge(Vertex *tail, Vertex *head, double *weight, WtNetwork *nwp);
int WtFindithNonedge (Vertex *tail, Vertex *head, Dyad i, WtNetwork *nwp);
int WtGetRandNonedge(Vertex *tail, Vertex *head, WtNetwork *nwp);
void Wtprintedge(Edge e, WtTreeNode *edges);
void WtInOrderTreeWalk(WtTreeNode *edges, Edge x);
void WtNetworkEdgeList(WtNetwork *nwp);
void WtShuffleEdges(Vertex *tails, Vertex *heads, double *weights, Edge nedges);
void WtDetShuffleEdges(Vertex *tails, Vertex *heads, double *weights, Edge nedges);
void WtDetUnShuffleEdges(Vertex *tails, Vertex *heads, double *weights, Edge nedges);

/* Others... */
Edge WtDesignMissing (Vertex a, Vertex b, WtNetwork *mnwp);
Edge WtEdgeTree2EdgeList(Vertex *tails, Vertex *heads, double *weights, WtNetwork *nwp, Edge nmax);

#endif
