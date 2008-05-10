/*
 *  File ergm/src/geodist.h
 *  Part of the statnet package, http://statnetproject.org
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) in
 *    http://statnetproject.org/attribution
 *
 * Copyright 2003 Mark S. Handcock, University of Washington
 *                David R. Hunter, Penn State University
 *                Carter T. Butts, University of California - Irvine
 *                Steven M. Goodreau, University of Washington
 *                Martina Morris, University of Washington
 * Copyright 2007 The statnet Development Team
 */

#ifndef GEODIST_H
#define GEODIST_H

#include "edgetree.h"

void node_geodesics (int *edgelist, int *nnodes, int *nodelist,
                     int *nedges, int *nodecolor, int *dist, 
                     int *Q, int *source);

void full_geodesic_distribution (int *edgelist, int *nnodes,
				 int *nodelist, int *nedges, 
				 int *nodecolor, int *dist, int *Q,
				 int *geodist);

void geodesic_matrix (int *edgelist, int *nnodes,
		      int *nodelist, int *nedges, 
		      int *nodecolor, int *distmat, int *Q);

void pair_geodesic (int *edgelist, int *nnodes, int *nodelist,
                     int *nedges, int *nodecolor, int *dist, 
		    int *Q, int *source, int *destination);

#endif
