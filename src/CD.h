/*  File src/CD.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2019 Statnet Commons
 */
#ifndef CD_H
#define CD_H

#include "MCMC.h"

/* *** don't forget tail-> head, so this function accepts tails first, not heads  */

SEXP CD_wrapper(// Network settings
                SEXP dn, SEXP dflag, SEXP bipartite,
                // Model settings
                SEXP nterms, SEXP funnames,
                SEXP sonames,
                // Proposal settings
                SEXP MHProposaltype, SEXP MHProposalpackage,
                SEXP attribs, SEXP maxout, SEXP maxin, SEXP minout,
                SEXP minin, SEXP condAllDegExact, SEXP attriblength,
                // Numeric inputs
                SEXP inputs,
                // Network state
                SEXP nedges,
                SEXP tails, SEXP heads,
                // MCMC settings
                SEXP eta, SEXP samplesize, 
                SEXP CDparams,
                SEXP verbose);
MCMCStatus CDSample(ErgmState *s,
		    double *eta, double *networkstatistics, 
		    int samplesize, int *CDparams,
                    Vertex *undotail, Vertex *undohead, double *extraworkspace, int verbose);
MCMCStatus CDStep(ErgmState *s,
		  double *eta, double *networkstatistics,
		  int *CDparams, int *staken,
		  Vertex *undotail, Vertex *undohead, double* extraworkspace,
		  int verbose);
#endif
