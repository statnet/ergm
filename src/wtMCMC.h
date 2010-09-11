#ifndef WTMCMC_H
#define WTMCMC_H

#include "wtedgetree.h"
#include "wtchangestat.h"
#include "wtMHproposal.h"
#include "wtmodel.h"

void WtMCMC_wrapper (int *heads, int *tails, double *weights, int *dnedges, double *baseline_weight,
                   int *maxpossibleedges,
                   int *dn, int *dflag, int *bipartite, 
                   int *nterms, char **funnames,
                   char **sonames, 
                   char **MHproposaltype, char **MHproposalpackage,
                   double *inputs, double *theta0, int *samplesize, 
                   double *sample, int *burnin, int *interval,  
                   int *newnetworkheads, 
                   int *newnetworktails, 
		     double *newnetworkweights,
                   int *fVerbose, 
                   int *maxedges);
void WtMCMCSample (char *MHproposaltype, char *MHproposalpackage,
		 double *theta, double *networkstatistics, 
		 long int samplesize, long int burnin, 
		 long int interval, int fVerbose,
	       	 WtNetwork *nwp, WtModel *m);
void WtMetropolisHastings (WtMHproposal *MHp,
			 double *theta, double *statistics, 
			 long int nsteps, long int *staken,
			 int fVerbose,
			 WtNetwork *nwp, WtModel *m);

#endif
