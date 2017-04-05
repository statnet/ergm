/*
 *  File ergm/src/geodist.h
 *  Part of the statnet package, http://statnet.org
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) in
 *    http://statnet.org/attribution
 *
 *  Copyright 2012 the statnet development team
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
