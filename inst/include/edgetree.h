/*  File inst/include/edgetree.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2018 Statnet Commons
 */
#ifndef EDGETREE_H
#define EDGETREE_H

#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#define MIN(a,b) ((a)<(b) ? (a) : (b))
#define MAX(a,b) ((a)<(b) ? (b) : (a))
#define DYADCOUNT(nnodes, bipartite, directed) ((bipartite)? (unsigned long)((nnodes)-(bipartite))*(unsigned long)(bipartite) : ((directed)? (unsigned long)(nnodes)*(unsigned long)((nnodes)-1) : (((unsigned long)(nnodes)*(unsigned long)((nnodes)-1))/2)))

/*typedef unsigned int Vertex; */
typedef int Vertex;
typedef unsigned int Edge;
typedef unsigned long int Dyad;

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

/* Dur_Inf is a structure containing information about durations of
edges in a network structure.
*/ 
typedef struct Dur_Infstruct {
  int time;
  int *lasttoggle;
} Dur_Inf;



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
Network NetworkInitialize(Vertex *tails, Vertex *heads, Edge nedges,
			  Vertex nnodes, int directed_flag, Vertex bipartite,
			  int lasttoggle_flag, int time, int *lasttoggle);
void NetworkDestroy(Network *nwp);
Network NetworkInitializeD(double *tails, double *heads, Edge nedges,
			   Vertex nnodes, int directed_flag, Vertex bipartite,
			   int lasttoggle_flag, int time, int *lasttoggle);

Network *NetworkCopy(Network *dest, Network *src);

/* Accessors. */
Edge EdgetreeSearch (Vertex a, Vertex b, TreeNode *edges);
Edge EdgetreeSuccessor (TreeNode *edges, Edge x);
Edge EdgetreePredecessor (TreeNode *edges, Edge x);
Edge EdgetreeMinimum (TreeNode *edges, Edge x);
Edge EdgetreeMaximum (TreeNode *edges, Edge x);

/* Modifiers. */

/* *** don't forget,  tails -> heads, so all the functions below using
   heads & tails, now list tails before heads */

int ToggleEdge (Vertex tail, Vertex head, Network *nwp);
int ToggleEdgeWithTimestamp (Vertex tail, Vertex head, Network *nwp);
int AddEdgeToTrees(Vertex tail, Vertex head, Network *nwp);
void AddHalfedgeToTree (Vertex a, Vertex b, TreeNode *edges, Edge *last_edge);
void CheckEdgetreeFull (Network *nwp);
int DeleteEdgeFromTrees(Vertex tail, Vertex head, Network *nwp);
int DeleteHalfedgeFromTree(Vertex a, Vertex b, TreeNode *edges,
		     Edge *last_edge);
void RelocateHalfedge(Edge from, Edge to, TreeNode *edges);

/* Duration functions. */
int ElapsedTime(Vertex tail, Vertex head, Network *nwp);
void TouchEdge(Vertex tail, Vertex head, Network *nwp);

/* Utility functions. */
int FindithEdge(Vertex *tail, Vertex *head, Edge i, Network *nwp);
int GetRandEdge(Vertex *tail, Vertex *head, Network *nwp);
// This one is implemented as a macro, since it's very simple and works exactly the same for weighted and unweighted.
#define GetRandDyad(tail, head, nwp)					\
  if((nwp)->bipartite){							\
    *(tail) = 1 + unif_rand() * (nwp)->bipartite;			\
    *(head) = 1 + (nwp)->bipartite + unif_rand() * ((nwp)->nnodes - (nwp)->bipartite); \
  }else{								\
    *(tail) = 1 + unif_rand() * (nwp)->nnodes;				\
    *(head) = 1 + unif_rand() * ((nwp)->nnodes-1);			\
    if(*(head)>=*(tail)) (*(head))++;					\
    									\
    if (!(nwp)->directed_flag && *(tail) > *(head)) {			\
      Vertex tmp = *(tail);						\
      *(tail) = *(head);						\
      *(head) = tmp;							\
    }									\
  }
int FindithNondge(Vertex *tail, Vertex *head, Dyad i, Network *nwp);
int GetRandNonedge(Vertex *tail, Vertex *head, Network *nwp);
void printedge(Edge e, TreeNode *edges);
void InOrderTreeWalk(TreeNode *edges, Edge x);
void NetworkEdgeList(Network *nwp);
void ShuffleEdges(Vertex *tails, Vertex *heads, Edge nedges);

/* Others... */
Edge DesignMissing (Vertex a, Vertex b, Network *mnwp);
Edge EdgeTree2EdgeList(Vertex *tails, Vertex *heads, Network *nwp, Edge nmax);

#endif
