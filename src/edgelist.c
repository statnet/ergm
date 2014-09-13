/*  File src/edgelist.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2013 Statnet Commons
 */
#include "edgelist.h"
/*********************
 unsigned int dEdgeListSearch

 A function to check whether a given edge is in an edgelist formatted
 as documented in edgelist.h using an array of double.

 Returns the index of the edge (counting from 1) if found, 0 if not.
*********************/

unsigned int dEdgeListSearch(Vertex tail, Vertex head, double *el){
  unsigned int nedges=el[0];
  unsigned int u=nedges,l=1;
  double *tails = el, *heads = el+nedges;

  if(nedges==0) return(0);

  while(l<u){
    unsigned int m = l + (u-l)/2;
    if(tail>tails[m] || (tail==tails[m] && head>heads[m])) l = m+1;
    else u = m;
  }

  if((u==l) && tail==tails[l] && head==heads[l]) return(l); else return(0);
}

/*********************
 unsigned int iEdgeListSearch

 A function to check whether a given edge is in an edgelist formatted
 as documented in edgelist.h using an array of int.

 Returns the index of the edge (counting from 1) if found, 0 if not.
*********************/

unsigned int iEdgeListSearch(Vertex tail, Vertex head, int *el){
  unsigned int nedges=el[0];
  unsigned int u=nedges,l=1;
  int *tails = el, *heads = el+nedges;

  if(nedges==0) return(0);

  while(l<u){
    unsigned int m = l + (u-l)/2;
    if(tail>tails[m] || (tail==tails[m] && head>heads[m])) l = m+1;
    else u = m;
  }

  if((u==l) && tail==tails[l] && head==heads[l]) return(l); else return(0);
}
