/*  File src/wtCD.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2019 Statnet Commons
 */
#ifndef WTCD_H
#define WTCD_H

#include "wtMCMC.h"

/* *** don't forget tail-> head, so this function accepts tails first, not heads  */
SEXP WtCD_wrapper(// Network settings
                  SEXP dn, SEXP dflag, SEXP bipartite,
                  // Model settings
                  SEXP nterms, SEXP funnames,
                  SEXP sonames,
                  // Proposal settings
                  SEXP MHProposaltype, SEXP MHProposalpackage,
                  // Numeric inputs
                  SEXP inputs,
                  // Network state
                  SEXP nedges,
                  SEXP tails, SEXP heads, SEXP weights,
                  // MCMC settings
                  SEXP eta, SEXP samplesize, 
                  SEXP CDparams,
                  SEXP verbose);
WtMCMCStatus WtCDSample(ErgmWtState *s,
                        double *eta, double *networkstatistics, 
			int samplesize, int *CDparams,
                        Vertex *undotail, Vertex *undohead, double *undoweight, double *extraworkspace,
                        int verbose);
WtMCMCStatus WtCDStep(ErgmWtState *s,
                      double *eta, double *networkstatistics,
                      int *CDparams, int *staken,
                      Vertex *undotail, Vertex *undohead, double *undoweight, double *extraworkspace,
                      int verbose);

#endif
