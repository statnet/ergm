/*  File inst/include/wtMCMC.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2017 Statnet Commons
 */
#ifndef WTMCMC_H
#define WTMCMC_H

#include "wtedgetree.h"
#include "wtchangestat.h"
#include "wtMHproposal.h"
#include "wtmodel.h"

// TODO: This might be worth moving into a common "constants.h".
typedef enum WtMCMCStatus_enum {
  WtMCMC_OK = 0,
  WtMCMC_TOO_MANY_EDGES = 1,
  WtMCMC_MH_FAILED = 2
} WtMCMCStatus;

/* *** don't forget tail-> head, so this function accepts tails first, not heads  */
void WtMCMC_wrapper(int *dnumnets, int *nedges, 
		    int *tails, int *heads, double *weights,
		    int *dn, int *dflag, int *bipartite, 
		    int *nterms, char **funnames,
		    char **sonames, 
		    char **MHproposaltype, char **MHproposalpackage,
		    double *inputs, double *theta0, int *samplesize, 
		    double *sample, int *burnin, int *interval,  
		    int *newnetworktails, 
		    int *newnetworkheads, 
		    double *newnetworkweights,
		    int *fVerbose, 
		    int *maxedges,
		    int *status);
WtMCMCStatus WtMCMCSample(WtMHproposal *MHp,
			   double *theta, double *networkstatistics, 
			   int samplesize, int burnin, 
			   int interval, int fVerbose, int nmax,
			   WtNetwork *nwp, WtModel *m);
WtMCMCStatus WtMetropolisHastings(WtMHproposal *MHp,
				   double *theta, double *statistics, 
				   int nsteps, int *staken,
				   int fVerbose,
				   WtNetwork *nwp, WtModel *m);

#endif
