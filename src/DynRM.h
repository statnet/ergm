#ifndef DYNRM_H
#define DYNRM_H

#include "MCMCDyn.h"
#include "MHproposals.h"

void MCMCDynPhase12(// Observed network.
		    int *tails, int *heads, int *n_edges,
        int *maxpossibleedges,
		    int *n_nodes, int *dflag, int *bipartite, 
		    // Ordering of formation and dissolution.
		    int *order_code, 
		    // Formation terms and proposals.
		    int *F_nterms, char **F_funnames, char **F_sonames, 
		    char **F_MHproposaltype, char **F_MHproposalpackage,
		    double *F_inputs, double *theta0, 
		    // Formation parameter fitting.
		    double *init_dev, double *gain,
		    int *phase1n_base, int *phase2n_base, int *phase2sub,
		    // Dissolution terms and proposals.
		    int *D_nterms, char **D_funnames, char **D_sonames, 
		    char **D_MHproposaltype, char **D_MHproposalpackage,
		    double *D_inputs, double *gamma0,
		    // Degree bounds.
		    int *attribs, int *maxout, int *maxin, int *minout,
		    int *minin, int *condAllDegExact, int *attriblength,
		    // MCMC settings.
		    int *burnin, int *interval, int *dyninterval,
		    // Space for output.
		    int *maxedges,
		    // Verbosity.
		    int *fVerbose);

void MCMCSampleDynPhase12(// Observed and discordant network.
			  Network *nwp,
			  // Ordering of formation and dissolution.
			  DynamOrder order,
			  // Formation terms and proposals.
			  Model *F_m, MHproposal *F_MH,
			  double *theta, 
			  // Formation parameter fitting.
			  double *dev, // DEViation of the current network's formation statistics from the target statistics.
			  double gain,
			  int phase1n_base, int phase2n_base, int phase2sub,
			  // Dissolution terms and proposals.
			  Model *D_m, MHproposal *D_MH,
			  double *gamma, 
			  // Dissolution parameter fitting --- to add later? -PK
			  // Space for output.
			  Edge nmax,
			  Vertex *difftime, Vertex *difftail, Vertex *diffhead,
			  // MCMC settings.
			  unsigned int burnin, unsigned int interval, unsigned int dyninterval,
			  // Verbosity.
			  int fVerbose);
#endif
