#ifndef MCMCDYN_H
#define MCMCDYN_H

#include "edgetree.h"
#include "MHproposal.h"
#include "changestat.h"
#include "model.h"

void MCMCDyn_init_common(int *tails, int *heads, int n_edges,
          int maxedges,
				  int n_nodes, int dflag, int bipartite, Network *nw,

				  int F_nterms, char *F_funnames, char *F_sonames, double *F_inputs, Model **F_m,
				  int D_nterms, char *D_funnames, char *D_sonames, double *D_inputs, Model **D_m,
				  
				  int *attribs, int *maxout, int *maxin, int *minout,
				  int *minin, int condAllDegExact, int attriblength, 

				  char *F_MHproposaltype, char *F_MHproposalpackage, MHproposal *F_MH,
				  char *D_MHproposaltype, char *D_MHproposalpackage, MHproposal *D_MH,
			 int fVerbose);


void MCMCDyn_finish_common(Network *nw,
			   Model *F_m,
			   Model *D_m,
			   MHproposal *F_MH,
			   MHproposal *D_MH);

void MCMCDyn_wrapper(// Starting network.
		     int *tails, int *heads, int *n_edges, int *maxpossibleedges,
		     int *dn, int *dflag, int *bipartite,
		     // Formation terms and proposals.
		     int *F_nterms, char **F_funnames, char **F_sonames, 
		     char **F_MHproposaltype, char **F_MHproposalpackage,
		     double *F_inputs, double *theta, 
		     // Dissolution terms and proposals.
		     int *D_nterms, char **D_funnames, char **D_sonames,
		     char **D_MHproposaltype, char **D_MHproposalpackage,
		     double *D_inputs, double *gamma0,
		     // Degree bounds.
		     int *attribs, int *maxout, int *maxin, int *minout,
		     int *minin, int *condAllDegExact, int *attriblength, 
		     // MCMC settings.
		     double *nsteps,  int *MH_interval,
		     double *burnin, double *interval,  
		     // Space for output.
		     double *F_sample, double *D_sample, 
		     int *newnetworktail, int *newnetworkhead, 
		     double *maxedges,
		     int *diffnetworktime, int *diffnetworktail, int *diffnetworkhead, 
		     // Verbosity.
		     int *fVerbose);

void MCMCSampleDyn(// Observed and discordant network.
		   Network *nwp,
		   // Formation terms and proposals.
		   Model *F_m, MHproposal *F_MH, double *theta,
		   // Dissolution terms and proposals.
		   Model *D_m, MHproposal *D_MH, double *gamma,
		   // Space for output.
		   double *F_stats, double *D_stats,// Do we still need these?
		   Edge nmax,
		   Vertex *difftime, Vertex *difftail, Vertex *diffhead,		    
		   // MCMC settings.
		   unsigned int nsteps, unsigned int MH_interval,
		   unsigned int burnin, unsigned int interval, 
		   // Verbosity.
		   int fVerbose);

void MCMCDyn1Step(Network *nwp,
		  Model *F_m, MHproposal *F_MH, double *theta,
		  Model *D_m, MHproposal *D_MH, double *gamma,
		  unsigned log_toggles,
		  double *F_stats, double *D_stats,
		  unsigned int nmax, Edge *nextdiffedge,
		  Vertex *difftime, Vertex *difftail, Vertex *diffhead,
		  unsigned int MH_interval,
		  int fVerbose);
#endif
