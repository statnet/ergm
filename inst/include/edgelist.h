/*  File inst/include/edgelist.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2018 Statnet Commons
 */
#ifndef EDGELIST_H
#define EDGELIST_H

#include "edgetree.h"

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


unsigned int dEdgeListSearch(Vertex tail, Vertex head, double *el);
unsigned int iEdgeListSearch(Vertex tail, Vertex head, int *el);

#endif
