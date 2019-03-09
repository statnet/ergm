/*  File src/wtnetstats.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2019 Statnet Commons
 */
#ifndef WTNETSTATS_H
#define WTNETSTATS_H

#include "ergm_wtedgetree.h"
#include "ergm_wtmodel.h"
#include "ergm_wtMHproposal.h"

void wt_network_stats_wrapper(int *tails, int *heads, double *weights, int *timings, int *time, int *lasttoggle, int *dnedges,
			   int *dn, int *dflag,  int *bipartite,
			   int *nterms, char **funnames,
			   char **sonames, double *inputs,  double *stats);
void WtSummStats(Edge n_edges, Vertex *tails, Vertex *heads, double *weights,
	       WtNetwork *nwp, WtModel *m, double *stats);
#endif
