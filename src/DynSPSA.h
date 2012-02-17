#ifndef DYNSPSA_H
#define DYNSPSA_H

#include "MCMCDyn.h"
#include "MHproposals.h"

void MCMCDynSPSA_wrapper(// Observed network.
		    int *tails, int *heads, int *n_edges,
		    int *n_nodes, int *dflag, int *bipartite, 
		    // Formation terms and proposals.
		    int *F_nterms, char **F_funnames, char **F_sonames, int *F_offset,
		    char **F_MHproposaltype, char **F_MHproposalpackage,
		    double *F_inputs, double *F_theta0, 
		    double *init_dev,
		    // SPSA settings.
		    double *a,
		    double *alpha,
		    double *A,
		    double *c,
		    double *gamma,
		    int *iterations,
		    // Dissolution terms and proposals.
		    int *D_nterms, char **D_funnames, char **D_sonames, 
		    char **D_MHproposaltype, char **D_MHproposalpackage,
		    double *D_inputs, double *D_theta0,
		    // Degree bounds.
		    int *attribs, int *maxout, int *maxin, int *minout,
		    int *minin, int *condAllDegExact, int *attriblength,
		    // MCMC settings.
		    int *burnin, int *interval, int *dyninterval,
		    // Space for output.
		    int *maxedges,
		    double *obj_history,
		    // Verbosity.
		    int *fVerbose);

double MCMCSampleDynObjective(Network *nwp,
			      // Formation terms and proposals.
			      Model *F_m, MHproposal *F_MH, double *F_theta,
			      // Dissolution terms and proposals.
			      Model *D_m, MHproposal *D_MH, double *D_theta,
			      double *F_stats, double *D_stats,
			      unsigned int nmax,
			      Vertex *difftime, Vertex *difftail, Vertex *diffhead,
			      // MCMC settings.
			      unsigned int dyninterval,
			      unsigned int burnin,
			      unsigned int S,
			      double *F_stats_acc, double* F_stats2_acc, int *use_var,
			      int fVerbose);

void MCMCDynSPSA(// Observed and discordant network.
		       Network *nwp,
		       // Formation terms and proposals.
		       Model *F_m, int *F_offset, MHproposal *F_MH,
		       double *F_theta, 
		       // Formation parameter fitting.
		       double *dev, // DEViation of the current network's formation statistics from the target statistics.
		       // SPSA settings.
		       double a,
		       double alpha,
		       double A,
		       double c,
		       double gamma,
		       unsigned int iterations,
		       // Dissolution terms and proposals.
		       Model *D_m, MHproposal *D_MH,
		       double *D_theta, 
		       // Dissolution parameter fitting --- to add later? -PK
		       // Space for output.
		       Edge nmax,
		       double *obj_history,
		       Vertex *difftime, Vertex *difftail, Vertex *diffhead,
		       // MCMC settings.
		       unsigned int burnin, unsigned int interval, unsigned int dyninterval,
		       // Verbosity.
		       int fVerbose);


#endif
