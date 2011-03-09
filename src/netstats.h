#ifndef NETSTATS_H
#define NETSTATS_H

#include "edgetree.h"
#include "model.h"
#include "MHproposal.h"

/* *** don't forget tail -> head, so these functions accept tails first, not heads */

void network_stats_wrapper(int *tails, int *heads, int *dnedges, 
			   int *dn, int *dflag,  int *bipartite,
			   int *nterms, char **funnames,
			   char **sonames, double *inputs,  double *stats);
void SummStats(Edge n_edges, Vertex *tails, Vertex *heads,
	       Network *nwp, Model *m, double *stats);
#endif
