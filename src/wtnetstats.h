#ifndef WTNETSTATS_H
#define WTNETSTATS_H

#include "wtedgetree.h"
#include "wtmodel.h"
#include "wtMHproposal.h"

void wt_network_stats_wrapper(int *tails, int *heads, double *weights, int *timings, int *time, int *lasttoggle, int *dnedges,
			   int *dn, int *dflag,  int *bipartite,
			   int *nterms, char **funnames,
			   char **sonames, double *inputs,  double *stats);
void WtSummStats(Edge n_edges, Vertex *tails, Vertex *heads, double *weights,
	       WtNetwork *nwp, WtModel *m, double *stats);
#endif
