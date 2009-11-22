#ifndef PILA_H
#define PILA_H

#include "edgetree.h"
#include "changestat.h"
#include "MHproposal.h"
#include "model.h"

void PILA_wrapper(int *heads, int *tails, int *dnedges,
		  int *dn, int *dflag, int *bipartite, 
                   int *nterms, char **funnames,
                   char **sonames, 
                   char **MHproposaltype, char **MHproposalpackage,
                   double *inputs, double *theta,
		   int *samplesize, int *interval, 
                   double *sample, int *burnin, 
		   double *alpha, double *gamma,
		  int *fVerbose, 
                   int *attribs, int *maxout, int *maxin, int *minout,
                   int *minin, int *condAllDegExact, int *attriblength, 
		  int *mheads, int *mtails, int *mdnedges,
		  double *theta_mean_save, double *XtX_save, double *XtY_save, double *beta_save, 
		  double *direction_save, double *dtheta_save,
		  int *insensitive_save, int *ineffectual_save, int *dropped_save);

void PILASample (char *MHproposaltype, char *MHproposalpackage,
		  double *theta, double *networkstatistics, 
		 long int samplesize, unsigned int interval,
		 long int burnin,
		 double alpha, double gamma,
		  int hammingterm, int fVerbose,
		 Network *nwp, Model *m, DegreeBound *bd,
		 double *theta_mean_save, double *XtX_save, double *XtY_save, double *beta_save, 
		 double *direction_save, double *dtheta_save,
		 int *insensitive_save, int *ineffectual_save, int *dropped_save);

#endif
