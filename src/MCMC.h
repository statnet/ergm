/*  File inst/include/MCMC.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2019 Statnet Commons
 */
#ifndef _MCMC_H_
#define _MCMC_H_

#include "ergm_constants.h"
#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_MHproposal.h"
#include "ergm_model.h"
#include "ergm_state.h"

MCMCStatus MCMCSample(ErgmState *s,
		      double *eta, double *networkstatistics, 
		      int samplesize, int burnin, 
		      int interval, int nmax, int verbose);
MCMCStatus MetropolisHastings(ErgmState *s,
			      double *eta, double *statistics, 
			      int nsteps, int *staken,
			      int verbose);

MCMCStatus MCMCSamplePhase12(ErgmState *s,
  double *eta, double gain, double *meanstats,
  int nphase1, int nsubphases, double *networkstatistics, 
  int samplesize, int burnin, 
  int interval, int verbose);

#endif
