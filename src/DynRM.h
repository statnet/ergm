#ifndef DYNRM_H
#define DYNRM_H

#include "MCMCDyn.h"
#include "MHproposals.h"

void MCMCDynRMPhase2_wrapper(// Observed network.
			     int *tails, int *heads, int *time, int *lasttoggle, int *n_edges,
			     int *n_nodes, int *dflag, int *bipartite, 
			     // Formation terms and proposals.
			     int *F_nterms, char **F_funnames, char **F_sonames,
			     char **F_MHproposaltype, char **F_MHproposalpackage,
			     double *F_inputs, double *F_eta0, 
			     // Dissolution terms and proposals.
			     int *D_nterms, char **D_funnames, char **D_sonames, 
			     char **D_MHproposaltype, char **D_MHproposalpackage,
			     double *D_inputs, double *D_eta0, 
			     // Parameter fitting.
			     int *M_nterms, char **M_funnames, char **M_sonames, double *M_inputs,
			     double *init_dev,
			     int *phase2n,
			     double *WinvGradient,
			     double *jitter, double *dejitter,
			     // Degree bounds.
			     int *attribs, int *maxout, int *maxin, int *minout,
			     int *minin, int *condAllDegExact, int *attriblength,
			     // MCMC settings.
			     int *RM_burnin, int *RM_interval, int *MH_interval,
			     // Space for output.
			     int *maxedges, int *maxchanges,
			     int *newnetworktail, int *newnetworkhead, 
			     double *opt_history,
			     // Verbosity.
			     int *fVerbose,
			     int *status);
MCMCDynStatus MCMCDynRMPhase2(// Observed and discordant network.
			      Network *nwp,
			      // Formation terms and proposals.
			      Model *F_m, MHproposal *F_MH,
			      double *F_eta, 
			      // Dissolution terms and proposals.
			      Model *D_m, MHproposal *D_MH,
			      double *D_eta, 
			      // Monitored terms.
			      Model *M_m,
			      double *dev, // DEViation of the current network's targeted statistics from the target statistics.
			      int phase2n,
			      double *WinvGradient, double *jitter, double *dejitter,
			      
			      // Space for output.
			      Edge maxedges, Edge maxchanges,
			      Vertex *difftime, Vertex *difftail, Vertex *diffhead,
			      double *opt_history,
			      // MCMC settings.
			      unsigned int RM_burnin, unsigned int RM_interval, unsigned int MH_interval,
			      // Verbosity.
			      int fVerbose);
#endif
