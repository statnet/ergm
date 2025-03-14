/*  File src/MCMC.h.template.do_not_include_directly.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */

#include "ergm_constants.h"

MCMCStatus EDGETYPE(MCMCSample)(EDGETYPE(ErgmState) *s,
			   double *eta, double *networkstatistics, 
			   int samplesize, int burnin, 
			   int interval, int nmax, int verbose);
MCMCStatus EDGETYPE(MetropolisHastings)(EDGETYPE(ErgmState) *s,
				   double *eta, double *statistics, 
				   int nsteps, int *staken,
				   int verbose);

MCMCStatus EDGETYPE(MCMCSamplePhase12)(EDGETYPE(ErgmState) *s,
                               double *eta, unsigned int n_param, double gain,
                               int nphase1, int nsubphases,
                               int min_iterations, int max_iterations,
                               int burnin,
                               int interval, int verbose);

