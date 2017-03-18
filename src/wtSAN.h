/*  File src/wtSAN.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2017 Statnet Commons
 */
#ifndef WTSAN_H
#define WTSAN_H

#include "wtedgetree.h"
#include "wtchangestat.h"
#include "wtMHproposal.h"
#include "wtmodel.h"
#include "wtMCMC.h"

void WtSAN_wrapper (int * dnumnets, int *nedges,
		    int *tails, int *heads, double *weights,
		    int *dn, int *dflag, int *bipartite, 
		    int *nterms, char **funnames,
		    char **sonames, 
		    char **MHproposaltype, char **MHproposalpackage,
		    double *inputs, double *theta0, double *tau, int *samplesize, 
		    double *sample, int *burnin, int *interval,  
		    int *newnetworktails, 
		    int *newnetworkheads, 
		    double *newnetworkweights,
		    double *invcov,
		    int *fVerbose, 
		    int *maxedges,
		    int *status);

WtMCMCStatus WtSANSample (WtMHproposal *MHp,
		double *theta, double *invcov, double *tau, double *networkstatistics, 
		int samplesize, int burnin, 
	        int interval, int fVerbose, int nmax,
		WtNetwork *nwp, WtModel *m);
WtMCMCStatus WtSANMetropolisHastings (WtMHproposal *MHp,
			 double *theta, double *invcov, double *tau, double *statistics, 
			 int nsteps, int *staken,
			 int fVerbose,
			 WtNetwork *nwp, WtModel *m);
#endif
