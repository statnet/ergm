/*  File src/CD.h.template.do_not_include_directly.h in package ergm, part of
 *  the Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
MCMCStatus ETYPE(CDSample)(ETYPE(ErgmState) *s,
                        double *eta, double *networkstatistics, 
			int samplesize, int *CDparams,
                        CD_UNDOS_RECEIVE, double *extraworkspace,
                        int verbose);
MCMCStatus ETYPE(CDStep)(ETYPE(ErgmState) *s,
                      double *eta, double *networkstatistics,
                      int *CDparams, int *staken,
                      CD_UNDOS_RECEIVE, double *extraworkspace,
                      int verbose);
