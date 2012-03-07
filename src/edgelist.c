/*
 *  File ergm/src/edgelist.c
 *  Part of the statnet package, http://statnet.org
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) in
 *    http://statnet.org/attribution
 *
 *  Copyright 2012 the statnet development team
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
  unsigned int u=nedges,l=1,m;

  if(nedges==0) return(0);

  do{
    m = l + (u-l)/2;
    if(tail==el[m] && head==el[nedges+m]) break;
    if(tail>el[m] || (tail==el[m] && head>el[nedges+m])) l = m+1;
    else u = m-1;
  }while(l<=u);

  if(l>u) return(0);
  else return(m);
}

/*********************
 unsigned int iEdgeListSearch

 A function to check whether a given edge is in an edgelist formatted
 as documented in edgelist.h using an array of int.

 Returns the index of the edge (counting from 1) if found, 0 if not.
*********************/

unsigned int iEdgeListSearch(Vertex tail, Vertex head, int *el){
  unsigned int nedges=el[0];
  unsigned int u=nedges,l=1,m;

  if(nedges==0) return(0);

  do{
    m = l + (u-l)/2;
    if(tail>el[m] || (tail==el[m] && head>el[nedges+m])) l = m+1;
    else u = m-1;
  }while(!(tail==el[m] && head==el[nedges+m]) && l<=u);

  if(l>u) return(0);
  else return(m);
}
