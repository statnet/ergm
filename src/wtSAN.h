#ifndef WTSAN_H
#define WTSAN_H

#include "wtedgetree.h"
#include "wtchangestat.h"
#include "wtMHproposal.h"
#include "wtmodel.h"

void WtSAN_wrapper (int * dnumnets, int *nedges,
		    int *heads, int *tails, double *weights,
		    int *maxpossibleedges,
		    int *dn, int *dflag, int *bipartite, 
		    int *nterms, char **funnames,
		    char **sonames, 
		    char **MHproposaltype, char **MHproposalpackage,
		    double *inputs, double *theta0, double *tau, int *samplesize, 
		    double *sample, int *burnin, int *interval,  
		    int *newnetworkheads, 
		    int *newnetworktails, 
		    double *newnetworkweights,
		    double *invcov,
		    int *fVerbose, 
		    int *maxedges);

void WtSANSample (char *MHproposaltype, char *MHproposalpackage,
		double *theta, double *invcov, double *tau, double *networkstatistics, 
		long int samplesize, long int burnin, 
		long int interval, int fVerbose,
		WtNetwork *nwp, WtModel *m);
void WtSANMetropolisHastings (WtMHproposal *MHp,
			 double *theta, double *invcov, double *tau, double *statistics, 
			 long int nsteps, long int *staken,
			 int fVerbose,
			 WtNetwork *nwp, WtModel *m);
#endif
