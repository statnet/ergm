/*  File src/wtCD.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2020 Statnet Commons
 */
#ifndef WTCD_H
#define WTCD_H

#include "wtMCMC.h"

/* *** don't forget tail-> head, so this function accepts tails first, not heads  */
void WtCD_wrapper(int *nedges, 
		    int *tails, int *heads, double *weights,
		    int *dn, int *dflag, int *bipartite, 
		    int *nterms, char **funnames,
		    char **sonames, 
		    char **MHProposaltype, char **MHProposalpackage,
		    double *inputs, double *theta0, int *samplesize, int *CDparams,
		    double *sample,
		    int *fVerbose, 
		    int *status);
MCMCStatus WtCDSample(WtMHProposal *MHp,
			double *theta, double *networkstatistics, 
			int samplesize, int *CDparams, Vertex *undotail, Vertex *undohead, double *undoweight,
			int fVerbose,
			WtNetwork *nwp, WtModel *m, double *extraworkspace);
MCMCStatus WtCDStep(WtMHProposal *MHp,
		      double *theta, double *statistics, 
		      int *CDparams, int *staken, Vertex *undotail, Vertex *undohead, double *undoweight,
		      int fVerbose,
		      WtNetwork *nwp, WtModel *m, double *extraworkspace);

#endif
