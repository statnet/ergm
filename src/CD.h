/*  File src/CD.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2020 Statnet Commons
 */
#ifndef CD_H
#define CD_H

#include "MCMC.h"

/* *** don't forget tail-> head, so this function accepts tails first, not heads  */

void CD_wrapper(int *dnedges,
		  int *tails, int *heads,
		  int *dn, int *dflag, int *bipartite, 
		  int *nterms, char **funnames,
		  char **sonames, 
		  char **MHProposaltype, char **MHProposalpackage,
		double *inputs, double *theta0, int *samplesize, int *CDparams,
		  double *sample, 
		  int *fVerbose, 
		  int *attribs, int *maxout, int *maxin, int *minout,
		  int *minin, int *condAllDegExact, int *attriblength, 
		int *status);
MCMCStatus CDSample(MHProposal *MHp,
		    double *theta, double *networkstatistics, 
		    int samplesize, int *CDparams, Vertex *undotails, Vertex *undoheads,
		    int fVerbose,
		    Network *nwp, Model *m, double *extraworkspace);
MCMCStatus CDStep(MHProposal *MHp,
		  double *theta, double *statistics,
		  int *CDparams, int *staken, Vertex *undotail, Vertex *undohead, 
		  int fVerbose,
		  Network *nwp, Model *m, double *extraworkspace);
#endif
