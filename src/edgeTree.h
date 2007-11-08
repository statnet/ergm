#ifndef EDGETREE_H
#define EDGETREE_H

#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#define MAXEDGES 100000
#define MIN(a,b) ((a)<(b) ? (a) : (b))
#define MAX(a,b) ((a)<(b) ? (b) : (a))

/*typedef unsigned int Vertex;
typedef unsigned int Edge; */
typedef int Vertex;
typedef int Edge;



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

/* WtTreeNode is just like TreeNode but with an extra field for a
   weight, or value, that might be associated with the node */
typedef struct WtTreeNodestruct {
  Vertex value;      /*  the vertex at the other end of the edge  */
  Edge parent;   /*  parent of this node in the tree (0 for root) */
  Edge left;     /*  left child (0 if none)  */
  Edge right;    /*  right child (0 if none) */
  double weight;
} WtTreeNode;

/* Dur_Inf is a structure containing information about durations of
edges in a network structure.
*/ 
typedef struct Dur_Infstruct {
  int MCMCtimer;
  int *lasttoggle;
/*  double mean_edge_duration; This is probably not a good idea */
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
   next_inedge and next_outedge are continually updated to give
     the smallest index of an edge object not being used.  
   outdegree[] and indegree[] are continually updated to give
     the appropriate degree values for each vertex.  These should
     point to Vertex-vectors of length nnodes.  
   value:  optional value(s) associated with this network 
*/
typedef struct Networkstruct {
  TreeNode *inedges;
  TreeNode *outedges;
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
} Network;



/* WtNetwork is just like Network except it is for a network with 
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
Network NetworkInitialize(Vertex *heads, Vertex *tails, Edge nedges,
			  Vertex nnodes, int directed_flag, Vertex bipartite,
			  int lasttoggle_flag);
Network NetworkInitializeD(double *heads, double *tails, Edge nedges,
			   Vertex nnodes, int directed_flag, Vertex bipartite,
			   int lasttoggle_flag);
WtNetwork WtNetworkInitialize(int *heads, int *tails, double *weights,
			      int nedges, int nnodes, int directed_flag,
			      int bipartite);

Edge DesignMissing (Vertex a, Vertex b, Network *mnwp);
Edge WtDesignMissing (Vertex a, Vertex b, WtNetwork *mnwp);

void NetworkDestroy(Network *nwp);
void WtNetworkDestroy(WtNetwork *nwp);

Edge EdgetreeSearch (Vertex a, Vertex b, TreeNode *edges);
Edge WtEdgetreeSearch (Vertex a, Vertex b, WtTreeNode *edges);

Edge EdgetreeSuccessor (TreeNode *edges, Edge x);
Edge WtEdgetreeSuccessor (WtTreeNode *edges, Edge x);

Edge EdgetreeMinimum (TreeNode *edges, Edge x);
Edge WtEdgetreeMinimum (WtTreeNode *edges, Edge x);

int ToggleEdge (Vertex head, Vertex tail, Network *nwp);
int WtToggleEdge (Vertex head, Vertex tail, double weight, WtNetwork *nwp);

int ToggleEdgeWithTimestamp (Vertex head, Vertex tail, Network *nwp);
int WtToggleEdgeWithTimestamp (Vertex head, Vertex tail, double weight, WtNetwork *nwp);

int ElapsedTime (Vertex head, Vertex tail, Network *nwp);

int AddEdgeToTrees(Vertex head, Vertex tail, Network *nwp);
int WtAddEdgeToTrees(Vertex head, Vertex tail, double weight, WtNetwork *nwp);

void AddHalfedgeToTree (Vertex a, Vertex b, TreeNode *edges, 
		   Edge *next_edge);
void WtAddHalfedgeToTree (Vertex a, Vertex b, double weight, 
			  WtTreeNode *edges, Edge *next_edge);

int DeleteEdgeFromTrees(Vertex head, Vertex tail, Network *nwp);
int WtDeleteEdgeFromTrees(Vertex head, Vertex tail, WtNetwork *nwp);

int DeleteHalfedgeFromTree(Vertex a, Vertex b, TreeNode *edges,
		     Edge *next_edge);
int WtDeleteHalfedgeFromTree(Vertex a, Vertex b, WtTreeNode *edges,
		     Edge *next_edge);

int FindithEdge (Vertex *head, Vertex *tail, Edge i, Network *nwp);
int WtFindithEdge (Vertex *head, Vertex *tail, Edge i, WtNetwork *nwp);

void printedge(Edge e, TreeNode *edges);
void Wtprintedge(Edge e, WtTreeNode *edges);

void InOrderTreeWalk(TreeNode *edges, Edge x);
void WtInOrderTreeWalk(WtTreeNode *edges, Edge x);

void NetworkEdgeList(Network *nwp);
void WtNetworkEdgeList(WtNetwork *nwp);

void TouchEdge(Vertex head, Vertex tail, Network *nwp);
Edge EdgeTree2EdgeList(Vertex *heads, Vertex *tails, Network *nwp, Edge nmax);

/* Below are some functions that only exist for weighted (valued) networks */

double EdgeWeight (Vertex head, Vertex tail, WtNetwork *nwp);

#endif









