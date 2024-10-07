/*  File src/geodist.c in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2024 Statnet Commons
 */
#include "geodist.h"

/* The geodist functions are based on breadth-first search.

From Cormen, Leiserson, and Rivest (1994):

 "To keep track of progress, breadth-first search colors each vertex white,
gray, or black.  All vertices start out white and may later become gray
and then black.  A vertex is DISCOVERED the first time it is encountered
during the search, at which time it becomes nonwhite.  Gray and black
vertices, therefore, have been discovered, but breadth-first search
distinguishes between them to ensure that the search proceeds in a
breadth-first manner.  If (u,v) is an edge and vertex u is black, then
vertex v is either gray or black; that is, all vertices adjacent to black
vertices have been discovered.  Gray vertices may have some adjacent white
vertices; they represent the frontier between discovered and undiscovered
vertices."

Note, however, that the distinction between gray and black is unimportant;
the queue always contains the gray nodes, in case we ever needed to know 
what they are.  Thus, the pseudocode below uses only two colors, white
and nonwhite.

We need 2 n-vectors:  color and distance.  We assume s is the source
vertex.  We also need a first-in, first-out queue of length at most n.

Initialize:
  set all color=white
  set all distance=infinity
  set color[s]=nonwhite
  set distance[s]=0
  set queue={s}

Main Loop:
  while queue is not empty  {
    pop the bottom value off the queue; call it u
    for each v adjacent to u {
      if color[v]=white {
        color[v]=nonwhite
        d[v]=1+d[u]
        push v on top of queue
      }
    }
  }
*/

/*********  FUNCTION node_geodesics
 All edges are assumed directed; thus, for undirected graphs, edgelist 
 must contain BOTH (i,j) and (j,i).  edgelist should be in the form of 
 a 2m-vector of type int, where the ith edge connects from node 
 edgelist[2*i-2] to node edgelist[2*i-1].  It is assumed that all of the
 edges originating from node u are listed consecutively beginning
 at edgelist [2*nodelist[u-1]]; thus, nodelist is a vector of length
 nnodes that must be prepared by the calling subroutine.  
 edgelist, nnodes, nodelist, nedges, and source are not altered by this 
 code; the remaining parameters (nodecolor, dist, Q) are.  The goal of 
 this subroutine is to correctly modify the dist vector for use by 
 other functions.  Upon return, the value of dist[i-1] should equal 
 the geodesic length from the source node to node i, or nnodes if no 
 path connects these two nodes.  */
void node_geodesics (int *edgelist, int *nnodes, int *nodelist,
                     int *nedges, int *nodecolor, int *dist, 
                     int *Q, int *source) {
  int i, j, u, v, n=*nnodes, twoe = 2*(*nedges), Qbottom=0, Qtop=0;

  for (i=0; i<n; i++) {
    nodecolor[i]=0; /* WHITE */
    dist[i]=n; /* Here, n means infinity */
  }
  nodecolor[*source-1]=1; /* NONWHITE */
  dist[*source-1]=0; 
  Q[Qtop++]=*source;  /* Push source onto top of queue */
  while (Qbottom<Qtop) {  /* Repeat until queue is empty */
    u=Q[Qbottom++]; /* Pop vertex off bottom of queue (it must be NONWHITE) */
    for (j=2*nodelist[u-1]; j<twoe && edgelist[j]==u; j+=2) {
      v=edgelist[j+1];
      if (nodecolor[v-1]==0) { /* WHITE */
        nodecolor[v-1]=1; /* NONWHITE */
        dist[v-1] = dist[u-1]+1; /* Node v is one step farther than node u */
        Q[Qtop++]=v;  /* Push v onto top of queue */
      }
    }
  }
}

/*************  FUNCTION full_geodesic_distribution
 See notes for node_geodesics about input parameters.  The additional
 input here is geodist, which is an n-vector whose values are modified
 by this function so that upon return:
       geodist[0]: # node pairs at a distance of infinity
       geodist[i], i>0:  # node pairs at a distance of i
 Simple check:  The sum of the entries of geodist should equal n(n-1). */
void full_geodesic_distribution (int *edgelist, int *nnodes,
				 int *nodelist, int *nedges, 
				 int *nodecolor, int *dist, int *Q,
				 int *geodist) {
  int i, j, n=*nnodes;

  for(i=0; i<n; i++)
    geodist[i]=0;
  /* Rprintf("nnodes = %d, nedges = %d \n", n, *nedges); */
  for(i=1; i<=n; i++) {
    node_geodesics(edgelist, nnodes, nodelist, nedges, nodecolor, dist,
		   Q, &i);
    for(j=0; j<i-1; j++) 
      ++geodist[dist[j]-1];
    for(j=i; j<n; j++) 
      ++geodist[dist[j]-1];
  }
}
