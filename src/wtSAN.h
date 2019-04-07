/*  File src/wtSAN.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2019 Statnet Commons
 */
#ifndef WTSAN_H
#define WTSAN_H

#include "ergm_wtedgetree.h"
#include "ergm_wtchangestat.h"
#include "ergm_wtMHproposal.h"
#include "ergm_wtmodel.h"
#include "wtMCMC.h"

void WtSAN_wrapper (int *nedges,
		    int *tails, int *heads, double *weights,
		    int *dn, int *dflag, int *bipartite, 
		    int *nterms, char **funnames,
		    char **sonames, 
		    char **MHProposaltype, char **MHProposalpackage,
		    double *inputs, double *tau,
		    double *sample, double *prop_sample,
		    int *samplesize, int *nsteps,  
		    int *newnetworktails, 
		    int *newnetworkheads, 
		    double *newnetworkweights,
		    double *invcov,
		    int *fVerbose, 
		    int *maxedges,
		    int *status);

WtMCMCStatus WtSANSample (WtMHProposal *MHp,
		double *invcov, double *tau, double *networkstatistics, double *prop_networkstatistics,
		int samplesize, int nsteps, 
	        int fVerbose, int nmax,
		WtNetwork *nwp, WtModel *m);
WtMCMCStatus WtSANMetropolisHastings (WtMHProposal *MHp,
			 double *invcov, double *tau, double *statistics, double *prop_statistics,
			 int nsteps, int *staken,
			 int fVerbose,
			 WtNetwork *nwp, WtModel *m);
#endif
