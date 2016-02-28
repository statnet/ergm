#ifndef CD_H
#define CD_H

#include "MCMC.h"

/* *** don't forget tail-> head, so this function accepts tails first, not heads  */

void CD_wrapper(int *dnumnets, int *dnedges,
		  int *tails, int *heads,
		  int *dn, int *dflag, int *bipartite, 
		  int *nterms, char **funnames,
		  char **sonames, 
		  char **MHproposaltype, char **MHproposalpackage,
		double *inputs, double *theta0, int *samplesize, int *CDparams,
		  double *sample, 
		  int *fVerbose, 
		  int *attribs, int *maxout, int *maxin, int *minout,
		  int *minin, int *condAllDegExact, int *attriblength, 
		int *status);
MCMCStatus CDSample(MHproposal *MHp,
		    double *theta, double *networkstatistics, 
		    int samplesize, int *CDparams, Vertex *undotails, Vertex *undoheads,
		    int fVerbose,
		    Network *nwp, Model *m, double *extraworkspace);
MCMCStatus CDStep(MHproposal *MHp,
		  double *theta, double *statistics,
		  int *CDparams, int *staken, Vertex *undotail, Vertex *undohead, 
		  int fVerbose,
		  Network *nwp, Model *m, double *extraworkspace);
#endif
