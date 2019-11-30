/*  File inst/include/wtMCMC.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2019 Statnet Commons
 */
#ifndef _ERGM_WTMCMC_H_
#define _ERGM_WTMCMC_H_

#include "ergm_wtedgetree.h"
#include "ergm_wtchangestat.h"
#include "ergm_wtMHproposal.h"
#include "ergm_wtmodel.h"
#include "ergm_wtstate.h"

// TODO: This might be worth moving into a common "constants.h".
typedef enum WtMCMCStatus_enum {
  WtMCMC_OK = 0,
  WtMCMC_TOO_MANY_EDGES = 1,
  WtMCMC_MH_FAILED = 2
} WtMCMCStatus;

/* *** don't forget tail-> head, so this function accepts tails first, not heads  */
WtMCMCStatus WtMCMCSample(ErgmWtState *s,
			   double *eta, double *networkstatistics, 
			   int samplesize, int burnin, 
			   int interval, int nmax, int verbose);
WtMCMCStatus WtMetropolisHastings(ErgmWtState *s,
				   double *eta, double *statistics, 
				   int nsteps, int *staken,
				   int verbose);

#endif
