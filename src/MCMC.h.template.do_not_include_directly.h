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

#define PROP_PRINT IFELSEEDGEWT(Rprintf("  (%d, %d) -> %f  ", MHp->toggletail[i], MHp->togglehead[i], MHp->toggleweight[i]), \
                                Rprintf("  (%d, %d)  ", MHp->toggletail[i], MHp->togglehead[i]))
#define PROP_CHANGESTATS EDGETYPE(ChangeStats)(MHp->ntoggles, MHp->toggletail, MHp->togglehead, IFEDGEWT(MHp->toggleweight,) nwp, m)
#define PROP_CHANGESTATS_DO EDGETYPE(ChangeStatsDo)(MHp->ntoggles, MHp->toggletail, MHp->togglehead, IFEDGEWT(MHp->toggleweight,) nwp, m)
#define PROP_CHANGESTATS_UNDO EDGETYPE(ChangeStatsUndo)(MHp->ntoggles, MHp->toggletail, MHp->togglehead, IFEDGEWT(MHp->toggleweight,) nwp, m)
#define PROP_COMMIT IFELSEEDGEWT(EDGETYPE(SetEdge)(MHp->toggletail[i], MHp->togglehead[i], MHp->toggleweight[i], nwp), \
                                 ToggleEdge(MHp->toggletail[i], MHp->togglehead[i], nwp))
#define PROP_FINISH IFELSEEDGEWT(EDGETYPE(SetEdge)(MHp->toggletail[MHp->ntoggles-1], MHp->togglehead[MHp->ntoggles-1], MHp->toggleweight[MHp->ntoggles-1], nwp), \
                                 ToggleEdge(MHp->toggletail[MHp->ntoggles-1], MHp->togglehead[MHp->ntoggles-1], nwp))
