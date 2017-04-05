#ifndef SAN_H
#define SAN_H

#include "edgetree.h"
#include "changestat.h"
#include "MHproposal.h"
#include "model.h"

/* *** don't forget tail -> head, so this function accepts tails first, not heads  */


void SAN_wrapper (int *dnumnets, int *nedges,
		  int *tails, int *heads,
		  int *maxpossibleedges,
		  int *dn, int *dflag, int *bipartite, 
		  int *nterms, char **funnames,
		  char **sonames, 
		  char **MHproposaltype, char **MHproposalpackage,
		  double *inputs, double *theta0, double *tau, int *samplesize, 
		  double *sample, int *burnin, int *interval,  
		  int *newnetworktails, 
		  int *newnetworkheads, 
		  double *invcov,
		  int *fVerbose, 
		  int *attribs, int *maxout, int *maxin, int *minout,
		  int *minin, int *condAllDegExact, int *attriblength, 
		  int *maxedges);

void SANSample (char *MHproposaltype, char *MHproposalpackage,
		double *theta, double *invcov, double *tau, double *networkstatistics, 
		long int samplesize, long int burnin, 
		long int interval, int hammingterm, int fVerbose,
		Network *nwp, Model *m, DegreeBound *bd);
void SANMetropolisHastings (MHproposal *MHp,
			 double *theta, double *invcov, double *tau, double *statistics, 
			 long int nsteps, long int *staken,
			 int hammingterm, int fVerbose,
			 Network *nwp, Model *m, DegreeBound *bd);
#endif
