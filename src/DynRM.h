#ifndef DYNRM_H
#define DYNRM_H

#include "MCMCDyn.h"
#include "MHproposals.h"

void MCMCDynRMPhase2_wrapper(// Observed network.
		    int *tails, int *heads, int *n_edges,
		    int *n_nodes, int *dflag, int *bipartite, 
		    // Formation terms and proposals.
		    int *F_nterms, char **F_funnames, char **F_sonames, int *F_offset,
		    char **F_MHproposaltype, char **F_MHproposalpackage,
		    double *F_inputs, double *F_theta0, 
		    // Formation parameter fitting.
		    double *init_dev,
		    int *phase2n,
		    // Dissolution terms and proposals.
		    int *D_nterms, char **D_funnames, char **D_sonames, 
		    char **D_MHproposaltype, char **D_MHproposalpackage,
		    double *D_inputs, double *D_theta0,
		    // Degree bounds.
		    int *attribs, int *maxout, int *maxin, int *minout,
		    int *minin, int *condAllDegExact, int *attriblength,
		    // MCMC settings.
		    int *RM_burnin, int *RM_interval, int *MH_interval,
		    double *invGradient,
		    // Space for output.
		    int *maxedges,
		    int *newnetworktail, int *newnetworkhead, 
		    double *obj_history,
		    // Verbosity.
		    int *fVerbose);
void MCMCDynRMPhase2(// Observed and discordant network.
			  Network *nwp,
			  // Formation terms and proposals.
			  Model *F_m, int *F_offset, MHproposal *F_MH,
			  double *F_theta, 
			  // Formation parameter fitting.
			  double *dev, // DEViation of the current network's formation statistics from the target statistics.
			  int phase2n,
			  // Dissolution terms and proposals.
			  Model *D_m, MHproposal *D_MH,
			  double *D_theta, 
			  // Dissolution parameter fitting --- to add later? -PK
			  // Space for output.
			  Edge nmax,
			  Vertex *difftime, Vertex *difftail, Vertex *diffhead,
			  double *obj_history,
			  // MCMC settings.
			  unsigned int RM_burnin, unsigned int RM_interval, unsigned int MH_interval, double *invGradient,
			  // Verbosity.
			  int fVerbose);
#endif
