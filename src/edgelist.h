#ifndef EDGELIST_H
#define EDGELIST_H

#include "edgetree.h"

/*
Edgelists are so simple, there is no point in definding a struct for
them. They are a serialization for a static network, with minimal
metadata, having the form

<type> edgelist[1+2*nedges] = {nedges,
                               head[1], head[2], ..., head[nedges],
                               tail[1], tail[2], ..., tail[nedges]}

Optionally, it may also include weights

<type> wtedgelist[1+3*nedges] = {nedges,
                                 head[1], head[2], ..., head[nedges],
                                 tail[1], tail[2], ..., tail[nedges],
                                 weight[1], weight[2], ..., weight[nedges]}

At this time, we shall require that the list shall be sorted by tails
(second column), with ties broken by heads. (Note that
as.matrix.network outputs edgelists sorted thus.)

*/


unsigned int dEdgeListSearch(Vertex head, Vertex tail, double *el);
unsigned int iEdgeListSearch(Vertex head, Vertex tail, int *el);

#endif
