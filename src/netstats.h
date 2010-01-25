#ifndef NETSTATS_H
#define NETSTATS_H

#include "edgetree.h"
#include "model.h"

void network_stats_wrapper(int *heads, int *tails, int *dnedges, 
			   int *dn, int *dflag,  int *bipartite,
			   int *nterms, char **funnames,
			   char **sonames, double *inputs,  double *stats);
void SummStats(Edge n_edges, Vertex *heads, Vertex *tails,
	       Network *nwp, Model *m, double *stats);
#endif
