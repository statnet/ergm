/*  File src/MPLEconddeg.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2014 Statnet Commons
 */
#ifndef MPLEconddeg_H
#define MPLEconddeg_H

#include "edgetree.h"
#include "changestat.h"
#include "MHproposal.h"
#include "model.h"

// TODO: This might be worth moving into a common "constants.h".
typedef enum MCMCStatus_enum {
  MCMC_OK = 0,
  MCMC_TOO_MANY_EDGES = 1,
  MCMC_MH_FAILED = 2
} MCMCStatus;

void MPLEconddeg_wrapper (int *dnumnets, int *dnedges,
		          int *tails, int *heads,
                   int *dn, int *dflag, int *bipartite, 
                   int *nterms, char **funnames,
                   char **sonames, 
                   char **MHproposaltype, char **MHproposalpackage,
                   double *inputs, double *theta0, int *samplesize, 
                   double *sample, int *burnin, int *interval,  
                   int *newnetworktails, 
                   int *newnetworkheads, 
                   int *fVerbose, 
                   int *attribs, int *maxout, int *maxin, int *minout,
                   int *minin, int *condAllDegExact, int *attriblength, 
                   int *maxedges,
		   int *status);
MCMCStatus CondDegSampler (MHproposal *MHp,
		 double *theta, double *networkstatistics, 
		 int samplesize, int burnin, 
		 int interval, int fVerbose, int nmax,
	       	 Network *nwp, Model *m);
#endif
