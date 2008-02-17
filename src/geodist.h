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
