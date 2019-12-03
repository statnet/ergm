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

MCMCStatus WtCDSample(WtErgmState *s,
                        double *eta, double *networkstatistics, 
			int samplesize, int *CDparams,
                        Vertex *undotail, Vertex *undohead, double *undoweight, double *extraworkspace,
                        int verbose);
MCMCStatus WtCDStep(WtErgmState *s,
                      double *eta, double *networkstatistics,
                      int *CDparams, int *staken,
                      Vertex *undotail, Vertex *undohead, double *undoweight, double *extraworkspace,
                      int verbose);

#endif
