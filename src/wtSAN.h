#ifndef WTSAN_H
#define WTSAN_H

#include "wtedgetree.h"
#include "wtchangestat.h"
#include "wtMHproposal.h"
#include "wtmodel.h"

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
		    int *maxedges);

void WtSANSample (WtMHproposal *MHp,
		double *theta, double *invcov, double *tau, double *networkstatistics, 
		int samplesize, int burnin, 
		int interval, int fVerbose,
		WtNetwork *nwp, WtModel *m);
void WtSANMetropolisHastings (WtMHproposal *MHp,
			 double *theta, double *invcov, double *tau, double *statistics, 
			 int nsteps, int *staken,
			 int fVerbose,
			 WtNetwork *nwp, WtModel *m);
#endif
