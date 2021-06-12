/*  File src/changestats.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2021 Statnet Commons
 */
#ifndef CHANGESTATS_H
#define CHANGESTATS_H

#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_Rutil.h"

Vertex CountTriangles (Vertex tail, Vertex head, int outcount,
                       int incount, Network *nwp);

void edgewise_path_recurse(Network *nwp, Vertex dest, Vertex curnode, 
     Vertex *visited, long int curlen, double *countv, long int maxlen, int semi);

void edgewise_cycle_census(Network *nwp, Vertex tail, Vertex head, 
                           double *countv, long int maxlen, int semi);

#endif
