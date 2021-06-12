/*  File inst/include/ergm_edgelist.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2021 Statnet Commons
 */
#ifndef _ERGM_EDGELIST_H_
#define _ERGM_EDGELIST_H_

#include "ergm_edgetree.h"

/*
Edgelists are so simple, there is no point in definding a struct for
them. They are a serialization for a static network, with minimal
metadata, having the form

<type> edgelist[1+2*nedges] = {nedges,
                               tail[1], tail[2], ..., tail[nedges],
                               head[1], head[2], ..., head[nedges]}

Optionally, it may also include weights:

<type> wtedgelist[1+3*nedges] = {nedges,
                                 tail[1], tail[2], ..., tail[nedges],
                                 head[1], head[2], ..., head[nedges],
                                 weight[1], weight[2], ..., weight[nedges]}

It doesn't affect the search functions. The list shall be sorted by
tails (first column), with ties broken by heads.

*/

/*********************
 unsigned int dEdgeListSearch

 A function to check whether a given edge is in an edgelist formatted
 as documented above, using an array of double.

 Returns the index of the edge (counting from 1) if found, 0 if not.
*********************/

static inline unsigned int dEdgeListSearch(Vertex tail, Vertex head, double *el){
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
 as documented above, using an array of int.

 Returns the index of the edge (counting from 1) if found, 0 if not.
*********************/

static inline unsigned int iEdgeListSearch(Vertex tail, Vertex head, int *el){
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

#endif
