#ifndef WTMCMC_H
#define WTMCMC_H

#include "wtedgetree.h"
#include "wtchangestat.h"
#include "wtMHproposal.h"
#include "wtmodel.h"

void WtMCMC_wrapper (int *dnumnets, int *nedges, 
		     int *tails, int *heads, double *weights,
		     int *maxpossibleedges,
		     int *dn, int *dflag, int *bipartite, 
		     int *nterms, char **funnames,
		     char **sonames, 
		     char **MHproposaltype, char **MHproposalpackage,
		     double *inputs, double *theta0, int *samplesize, 
		     double *sample, int *burnin, int *interval,  
		     int *newnetworktails, 
		     int *newnetworkheads, 
		     double *newnetworkweights,
		     int *fVerbose, 
		     int *maxedges);
void WtMCMCSample (WtMHproposal *MHp,
		 double *theta, double *networkstatistics, 
		 int samplesize, int burnin, 
		 int interval, int fVerbose,
	       	 WtNetwork *nwp, WtModel *m);
void WtMetropolisHastings (WtMHproposal *MHp,
			 double *theta, double *statistics, 
			 int nsteps, int *staken,
			 int fVerbose,
			 WtNetwork *nwp, WtModel *m);

#endif
