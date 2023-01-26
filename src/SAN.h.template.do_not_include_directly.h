/*  File src/SAN.h.template.do_not_include_directly.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2023 Statnet Commons
 */

MCMCStatus DISPATCH_SANSample(DISPATCH_ErgmState *s,
                       double *invcov, double *tau, double *networkstatistics, double *prop_networkstatistics,
                       int samplesize, int nsteps,
                       int nstats,
                       int *statindices,
                       int noffsets,
                       int *offsetindices,
                       double *offsets,
                       int verbose);
MCMCStatus DISPATCH_SANMetropolisHastings(DISPATCH_ErgmState *s,
                                   double *invcov, double *tau, double *statistics, double *prop_statistics,
                                   int nsteps, int *staken,
                                   int nstats,
                                   int *statindices,
                                   int noffsets,
                                   int *offsetindices,
                                   double *offsets,
                                          double *deltainvsig,
                                   int verbose);
