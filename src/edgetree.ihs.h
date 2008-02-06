#ifndef EDGETREE_IHS_H
#define EDGETREE_IHS_H

#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "edgeTree.h"



/* WtTreeNode is just like TreeNode but with an extra field for a
   weight, or value, that might be associated with the node */
typedef struct WtTreeNodestruct {
  Vertex value;      /*  the vertex at the other end of the edge  */
  Edge parent;   /*  parent of this node in the tree (0 for root) */
  Edge left;     /*  left child (0 if none)  */
  Edge right;    /*  right child (0 if none) */
  double weight;
} WtTreeNode;





/* Network is a structure containing all essential elements
   of a given network; it is a slightly rewritten version of the old Gptr,
   with some changes of awkard things, deletion of unnecessary things, and
   a new name more reflective of what it does!

   Some of the fields in a Network structure are:
   inedges and outedges are arrays of TreeNode that are used to 
     store all of the incoming and outgoing edges, respectively. 
   directed_flag is 1 or 0, depending on whether or not the 
     network is directed. 
   next_inedge and next_outedge are continually updated to give
     the smallest index of an edge object not being used.  
   outdegree[] and indegree[] are continually updated to give
     the appropriate degree values for each vertex.  These should
     point to Vertex-vectors of length nnodes.  
   value:  optional value(s) associated with this network 
   Dur_Inf:  See typedef for Dur_Infstruct above
*/
/* WtNetwork is a structure just like Network except it is for a network with 
   weighted (valued) edges.  */
typedef struct WtNetworkstruct {
  WtTreeNode *inedges;
  WtTreeNode *outedges;
  int directed_flag;
  Vertex bipartite; 
  Vertex nnodes;
  Edge nedges;
  Edge next_inedge;
  Edge next_outedge;
  Vertex *indegree;
  Vertex *outdegree;
  double *value;  
  Dur_Inf duration_info;  
} WtNetwork;

/* Note that there are really two versions of each function, one for 
   edges with weights (values) and one for those without.  This naturally
   results in a lot of redundant code.  In the future, if there are to be
   more types of edges/networks, it is probably worthwhile to write 
   single, more flexible versions of each function that are capable of 
   handling more than one type of network.  For now, this redundancy 
   approach, while it lacks elegance, has the advantage of speed. */
WtNetwork WtNetworkInitialize(int *heads, int *tails, double *weights,
			      int nedges, int nnodes, int directed_flag,
			      int bipartite);
Edge WtDesignMissing (Vertex a, Vertex b, WtNetwork *mnwp);
void WtNetworkDestroy(WtNetwork *nwp);
Edge WtEdgetreeSearch (Vertex a, Vertex b, WtTreeNode *edges);
Edge WtEdgetreeSuccessor (WtTreeNode *edges, Edge x);
Edge WtEdgetreeMinimum (WtTreeNode *edges, Edge x);
int WtToggleEdge (Vertex head, Vertex tail, double weight, WtNetwork *nwp);
int WtToggleEdgeWithTimestamp (Vertex head, Vertex tail, double weight, WtNetwork *nwp);
int ElapsedTime (Vertex head, Vertex tail, Network *nwp);
void TouchEdge(Vertex head, Vertex tail, Network *nwp);
int WtAddEdgeToTrees(Vertex head, Vertex tail, double weight, WtNetwork *nwp);
void WtAddHalfedgeToTree (Vertex a, Vertex b, double weight, 
			  WtTreeNode *edges, Edge *next_edge);
int WtDeleteEdgeFromTrees(Vertex head, Vertex tail, WtNetwork *nwp);
int WtDeleteHalfedgeFromTree(Vertex a, Vertex b, WtTreeNode *edges,
		     Edge *next_edge);
int WtFindithEdge (Vertex *head, Vertex *tail, Edge i, WtNetwork *nwp);
void Wtprintedge(Edge e, WtTreeNode *edges);
void WtInOrderTreeWalk(WtTreeNode *edges, Edge x);
void WtNetworkEdgeList(WtNetwork *nwp);

/* Below are some functions that only exist for weighted (valued) networks */

double EdgeWeight (Vertex head, Vertex tail, WtNetwork *nwp);

#endif









