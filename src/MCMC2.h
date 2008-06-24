#ifndef PILA_H
#define PILA_H

#include "edgetree.h"
#include "changestats.h"
#include "model.h"
#include "MCMC.h"


void PILA_wrapper(int *heads, int *tails, int *dnedges,
                   int *dn, int *dflag, int *bipartite, 
                   int *nterms, char **funnames,
                   char **sonames, 
                   char **MHproposaltype, char **MHproposalpackage,
                   double *inputs, double *theta, int *samplesize, 
                   int *interval, 
		   double *sample, int *burnin,
		   double *theta_burnin, double *sample_burnin,
		   double *alpha, double *gamma,
                   int *fVerbose, 
                   int *attribs, int *maxout, int *maxin, int *minout,
                   int *minin, int *condAllDegExact, int *attriblength, 
                   int *mheads, int *mtails, int *mdnedges);
void PILASample (char *MHproposaltype, char *MHproposalpackage,
		  double *theta, double *networkstatistics, 
		  long int samplesize, unsigned int interval, long int burnin,
		  double *theta_burnin, double *sample_burnin,
		  double alpha, double gamma,
		  int hammingterm, int fVerbose,
		  Network *nwp, Model *m, DegreeBound *bd);

#endif
