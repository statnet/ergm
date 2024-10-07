/*  File src/geodist.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2024 Statnet Commons
 */
#ifndef GEODIST_H
#define GEODIST_H

#include "ergm_edgetree.h"

void node_geodesics (int *edgelist, int *nnodes, int *nodelist,
                     int *nedges, int *nodecolor, int *dist, 
                     int *Q, int *source);

void full_geodesic_distribution (int *edgelist, int *nnodes,
				 int *nodelist, int *nedges, 
				 int *nodecolor, int *dist, int *Q,
				 int *geodist);

#endif
