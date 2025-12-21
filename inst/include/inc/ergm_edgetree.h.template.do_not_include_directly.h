/*  File inst/include/inc/ergm_edgetree.h.template.do_not_include_directly.h in
 *  package ergm, part of the Statnet suite of packages for network analysis,
 *  https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */

#include "ergm_edgetree_common.do_not_include_directly.h"

#ifdef __cplusplus
extern "C" {
#endif

/*  ETYPE(TreeNode) is a binary tree structure, which is how the edgelists
    are stored.  The root of the tree for vertex i will be inedges[i]
    or outedges[i].  inedges[0] and outedges[0] are unused, since the
    index 0 will indicate no link.  Indices are long unsigned integers,
    which means networks can contain 2^32-1= 4294967295 edges (enough to
    hold all edges in a 92682-vertex undirected network or a 65536-vertex
    directed network, assuming no multiple edges or loops).
*/
typedef struct ETYPE(TreeNodestruct) {
  Vertex value;      /*  the vertex at the other end of the edge  */
  Edge parent;   /*  parent of this node in the tree (0 for root) */
  Edge left;     /*  left child (0 if none)  */
  Edge right;    /*  right child (0 if none) */
  IFEWT(EWTTYPE weight;)
} ETYPE(TreeNode);

/* ETYPE(Network) is a structure containing all essential elements
   of a given network; it is a slightly rewritten version of the old Gptr,
   with some changes of awkard things, deletion of unnecessary things, and
   a new name more reflective of what it does!

   Some of the fields in a ETYPE(Network) structure are:
   inedges and outedges are arrays of TreeNode that are used to
     store all of the incoming and outgoing edges, respectively.
   directed_flag is 1 or 0, depending on whether or not the
     network is directed.
   last_inedge and last_outedge are continually updated to give
     the highest index of an edge object being used.
   outdegree[] and indegree[] are continually updated to give
     the appropriate degree values for each vertex.  These should
     point to Vertex-vectors of length nnodes+1.
*/
typedef struct ETYPE(Networkstruct) {
  ETYPE(TreeNode) *inedges;
  ETYPE(TreeNode) *outedges;
  Rboolean directed_flag;
  Vertex bipartite;
  Vertex nnodes;
  Edge nedges;
  Edge last_inedge;
  Edge last_outedge;
  Vertex *indegree;
  Vertex *outdegree;
  IFEWT(const char *eattrname;)
  Edge maxedges;

  unsigned int n_on_edge_change;
  unsigned int max_on_edge_change;
  void (**on_edge_change)(Vertex, Vertex, IFEWT(EWTTYPE,) void*, struct ETYPE(Networkstruct)*, EWTTYPE);
  void **on_edge_change_payload;
} ETYPE(Network);
typedef void (*ETYPE(On,NetworkEdgeChange))(Vertex, Vertex, IFEWT(EWTTYPE,) void*, ETYPE(Network)*, EWTTYPE);

/* Initialization and destruction. */

ETYPE(Network) *ETYPE(NetworkInitialize_noLT)(Vertex *tails, Vertex *heads, IFEWT(EWTTYPE *weights,) Edge nedges,
			       Vertex nnodes, Rboolean directed_flag, Vertex bipartite);
void ETYPE(NetworkDestroy)(ETYPE(Network) *nwp);

ETYPE(Network) *ETYPE(NetworkCopy)(ETYPE(Network) *src);

SEXP ETYPE(Network2Redgelist)(ETYPE(Network) *nwp);
ETYPE(Network) *ETYPE(Redgelist2,Network)(SEXP elR, Rboolean empty);

/* /\* Accessors. *\/ */
/* static inline Edge ETYPE(EdgetreeSearch) (Vertex a, Vertex b, ETYPE(TreeNode) *edges); */
/* static inline EWTTYPE ETYPE(GetEdge) (Vertex tail, Vertex head, ETYPE(Network) *nwp); */
/* static inline Edge ETYPE(EdgetreeSuccessor) (ETYPE(TreeNode) *edges, Edge x); */
/* static inline Edge ETYPE(EdgetreePredecessor) (ETYPE(TreeNode) *edges, Edge x); */
/* static inline Edge ETYPE(EdgetreeMinimum) (ETYPE(TreeNode) *edges, Edge x); */
/* static inline Edge ETYPE(EdgetreeMaximum) (ETYPE(TreeNode) *edges, Edge x); */

/* Modifiers. */

/* *** don't forget,  tails -> heads, so all the functions below using
   heads & tails, now list tails before heads */

void ETYPE(SetEdge) (Vertex tail, Vertex head, EWTTYPE weight, ETYPE(Network) *nwp);
int ETYPE(ToggleEdge) (Vertex tail, Vertex head, IFEWT(EWTTYPE weight,) ETYPE(Network) *nwp);
void ETYPE(AddEdgeToTrees)(Vertex tail, Vertex head, IFEWT(EWTTYPE weight,) ETYPE(Network) *nwp);
/* void ETYPE(AddHalfedgeToTree) (Vertex a, Vertex b, IFEWT(EWTTYPE weight,) ETYPE(TreeNode) *edges, Edge *last_edge); */
/* void ETYPE(CheckEdgetreeFull) (ETYPE(Network) *nwp); */
int ETYPE(DeleteEdgeFromTrees)(Vertex tail, Vertex head, ETYPE(Network) *nwp);
/* int ETYPE(DeleteHalfedgeFromTree)(Vertex a, Vertex b, ETYPE(TreeNode) *edges, */
/* 		     Edge *last_edge); */
/* void ETYPE(RelocateHalfedge)(Edge from, Edge to, ETYPE(TreeNode) *edges); */

/* Callback management. */
void ETYPE(AddOn,NetworkEdgeChange)(ETYPE(Network) *nwp, ETYPE(On,NetworkEdgeChange) callback, void *payload, unsigned int pos);
void ETYPE(DeleteOn,NetworkEdgeChange)(ETYPE(Network) *nwp, ETYPE(On,NetworkEdgeChange) callback, void *payload);

#ifdef __cplusplus
}
#endif

#include "ergm_edgetree_inline_template.do_not_include_directly.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Utility functions. */
int ETYPE(FindithEdge) (Vertex *tail, Vertex *head, IFEWT(EWTTYPE *weight,) Edge i, ETYPE(Network) *nwp);
int ETYPE(GetRandEdge)(Vertex *tail, Vertex *head, IFEWT(EWTTYPE *weight,) ETYPE(Network) *nwp);
int ETYPE(FindithNonedge) (Vertex *tail, Vertex *head, Dyad i, ETYPE(Network) *nwp);
int ETYPE(GetRandNonedge)(Vertex *tail, Vertex *head, ETYPE(Network) *nwp);
void ETYPE(InOrderTreeWalk)(ETYPE(TreeNode) *edges, Edge x);
void ETYPE(NetworkEdgeList)(ETYPE(Network) *nwp);
void ETYPE(ShuffleEdges)(Vertex *tails, Vertex *heads, IFEWT(EWTTYPE *weights,) Edge nedges);
void ETYPE(DetShuffleEdges)(Vertex *tails, Vertex *heads, IFEWT(EWTTYPE *weights,) Edge nedges);
void ETYPE(DetUnShuffleEdges)(Vertex *tails, Vertex *heads, IFEWT(EWTTYPE *weights,) Edge nedges);

/* Others... */
Edge ETYPE(EdgeTree2EdgeList)(Vertex *tails, Vertex *heads, IFEWT(EWTTYPE *weights,) ETYPE(Network) *nwp, Edge nmax);

#ifdef __cplusplus
}
#endif
