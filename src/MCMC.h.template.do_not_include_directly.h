/*  File src/MCMC.h.template.do_not_include_directly.h in package ergm, part of
 *  the Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */

#include "ergm_constants.h"

MCMCStatus ETYPE(MCMCSample)(ETYPE(ErgmState) *s,
			   double *eta, double *networkstatistics, 
			   int samplesize, int burnin, 
			   int interval, int nmax, int verbose);
MCMCStatus ETYPE(MetropolisHastings)(ETYPE(ErgmState) *s,
				   double *eta, double *statistics, 
				   int nsteps, int *staken,
				   int verbose);

MCMCStatus ETYPE(MCMCSamplePhase12)(ETYPE(ErgmState) *s,
                               double *eta, unsigned int n_param, double gain,
                               int nphase1, int nsubphases,
                               int min_iterations, int max_iterations,
                               int burnin,
                               int interval, int verbose);

#define PROP_PRINT IFELSEEWT(Rprintf("  (%d, %d) -> %f  ", MHp->toggletail[i], MHp->togglehead[i], MHp->toggleweight[i]), \
                                Rprintf("  (%d, %d)  ", MHp->toggletail[i], MHp->togglehead[i]))
#define PROP_CHANGESTATS ETYPE(ChangeStats)(MHp->ntoggles, MHp->toggletail, MHp->togglehead, IFEWT(MHp->toggleweight,) nwp, m)
#define PROP_CHANGESTATS_DO ETYPE(ChangeStatsDo)(MHp->ntoggles, MHp->toggletail, MHp->togglehead, IFEWT(MHp->toggleweight,) nwp, m)
#define PROP_CHANGESTATS_UNDO ETYPE(ChangeStatsUndo)(MHp->ntoggles, MHp->toggletail, MHp->togglehead, IFEWT(MHp->toggleweight,) nwp, m)
#define PROP_COMMIT IFELSEEWT(ETYPE(SetEdge)(MHp->toggletail[i], MHp->togglehead[i], MHp->toggleweight[i], nwp), \
                                 ToggleEdge(MHp->toggletail[i], MHp->togglehead[i], nwp))
#define PROP_FINISH IFELSEEWT(ETYPE(SetEdge)(MHp->toggletail[MHp->ntoggles-1], MHp->togglehead[MHp->ntoggles-1], MHp->toggleweight[MHp->ntoggles-1], nwp), \
                                 ToggleEdge(MHp->toggletail[MHp->ntoggles-1], MHp->togglehead[MHp->ntoggles-1], nwp))
