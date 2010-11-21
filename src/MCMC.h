#ifndef MCMC_H
#define MCMC_H

#include "edgetree.h"
#include "changestat.h"
#include "MHproposal.h"
#include "model.h"

void MCMC_wrapper (int *dnumnets, int *dnedges,
		   int *heads, int *tails,
                   int *maxpossibleedges,
                   int *dn, int *dflag, int *bipartite, 
                   int *nterms, char **funnames,
                   char **sonames, 
                   char **MHproposaltype, char **MHproposalpackage,
                   double *inputs, double *theta0, int *samplesize, 
                   double *sample, int *burnin, int *interval,  
                   int *newnetworkheads, 
                   int *newnetworktails, 
                   int *fVerbose, 
                   int *attribs, int *maxout, int *maxin, int *minout,
                   int *minin, int *condAllDegExact, int *attriblength, 
                   int *maxedges);
void MCMCSample (char *MHproposaltype, char *MHproposalpackage,
		 double *theta, double *networkstatistics, 
		 long int samplesize, long int burnin, 
		 long int interval, int hammingterm, int fVerbose,
	       	 Network *nwp, Model *m, DegreeBound *bd);
void MetropolisHastings (MHproposal *MHp,
			 double *theta, double *statistics, 
			 long int nsteps, long int *staken,
       int hammingterm,
       int fVerbose,
			 Network *nwp, Model *m, DegreeBound *bd);
void MCMCPhase12 (int *heads, int *tails, int *dnedges,
      int *maxpossibleedges,
		  int *dn, int *dflag, int *bipartite, 
		  int *nterms, char **funnames,
		  char **sonames, 
		  char **MHproposaltype, char **MHproposalpackage,
		  double *inputs, 
		  double *theta0, int *samplesize,
		  double *gain, double *meanstats, int *phase1, int *nsub,
		  double *sample, int *burnin, int *interval,  
		  int *newnetworkheads, 
		  int *newnetworktails, 
		  int *fVerbose, 
		  int *attribs, int *maxout, int *maxin, int *minout,
		  int *minin, int *condAllDegExact, int *attriblength, 
		  int *maxedges,
		  int *mheads, int *mtails, int *mdnedges);

void MCMCSamplePhase12 (char *MHproposaltype, char *MHproposalpackage,
  double *theta, double gain, double *meanstats,
  int nphase1, int nsubphases, double *networkstatistics, 
  long int samplesize, long int burnin, 
  long int interval, int hammingterm, int fVerbose,
  Network *nwp, Model *m, DegreeBound *bd);

#endif
