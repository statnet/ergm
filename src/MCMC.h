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

#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_MHproposal.h"
#include "ergm_model.h"
#include "ergm_state.h"

// TODO: This might be worth moving into a common "constants.h".
typedef enum MCMCStatus_enum {
  MCMC_OK = 0,
  MCMC_TOO_MANY_EDGES = 1,
  MCMC_MH_FAILED = 2
} MCMCStatus;

/* *** don't forget tail-> head, so this function accepts tails first, not heads  */

MCMCStatus MCMCSample(ErgmState *s,
		      double *theta, double *networkstatistics, 
		      int samplesize, int burnin, 
		      int interval, int fVerbose, int nmax);
MCMCStatus MetropolisHastings(ErgmState *s,
			      double *theta, double *statistics, 
			      int nsteps, int *staken,
			      int fVerbose);
void MCMCPhase12 (int *tails, int *heads, int *dnedges,
		  int *dn, int *dflag, int *bipartite, 
		  int *nterms, char **funnames,
		  char **sonames, 
		  char **MHProposaltype, char **MHProposalpackage,
		  double *inputs, 
		  double *theta0, int *samplesize,
		  double *gain, double *meanstats, int *phase1, int *nsub,
		  double *sample, int *burnin, int *interval,  
		  int *newnetworktails, 
		  int *newnetworkheads, 
		  int *fVerbose, 
		  int *attribs, int *maxout, int *maxin, int *minout,
		  int *minin, int *condAllDegExact, int *attriblength, 
		  int *maxedges,
		  int *mtails, int *mheads, int *mdnedges);

void MCMCSamplePhase12 (ErgmState *s,
  double *theta, double gain, double *meanstats,
  int nphase1, int nsubphases, double *networkstatistics, 
  int samplesize, int burnin, 
  int interval, int fVerbose);

#endif
