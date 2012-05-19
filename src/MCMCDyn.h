/*
 *  File ergm/src/MCMCDyn.h
 *  Part of the statnet package, http://statnet.org
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) in
 *    http://statnet.org/attribution
 *
 *  Copyright 2012 the statnet development team
 */
#ifndef MCMCDYN_H
#define MCMCDYN_H
#include<string.h>
#include "edgetree.h"
#include "MHproposal.h"
#include "changestat.h"
#include "model.h"

// TODO: This might be worth moving into a common "constants.h".
typedef enum MCMCDynStatus_enum {
  MCMCDyn_OK = 0,
  MCMCDyn_TOO_MANY_EDGES = 1,
  MCMCDyn_MH_FAILED = 2,
  MCMCDyn_TOO_MANY_CHANGES = 3
} MCMCDynStatus;


void MCMCDyn_init_common(int *tails, int *heads, int time, int *lasttoggle, int n_edges,
			 int n_nodes, int dflag, int bipartite, Network *nw,
			 
			 int F_nterms, char *F_funnames, char *F_sonames, double *F_inputs, Model **F_m,
			 int D_nterms, char *D_funnames, char *D_sonames, double *D_inputs, Model **D_m,
			 int M_nterms, char *M_funnames, char *M_sonames, double *M_inputs, Model **M_m,
			 
			 int *attribs, int *maxout, int *maxin, int *minout,
			 int *minin, int condAllDegExact, int attriblength,
			 
			 char *F_MHproposaltype, char *F_MHproposalpackage, MHproposal *F_MH,
			 char *D_MHproposaltype, char *D_MHproposalpackage, MHproposal *D_MH,
			 int fVerbose);


void MCMCDyn_finish_common(Network *nw,
			   Model *F_m,
			   Model *D_m,
			   Model *M_m,
			   MHproposal *F_MH,
			   MHproposal *D_MH);

void MCMCDyn_wrapper(// Starting network.
		     int *tails, int *heads, int *time, int *lasttoggle, int *n_edges,
		     int *n_nodes, int *dflag, int *bipartite,
		     // Formation terms and proposals.
		     int *F_nterms, char **F_funnames, char **F_sonames, 
		     char **F_MHproposaltype, char **F_MHproposalpackage,
		     double *F_inputs, double *F_eta, 
		     // Dissolution terms and proposals.
		     int *D_nterms, char **D_funnames, char **D_sonames,
		     char **D_MHproposaltype, char **D_MHproposalpackage,
		     double *D_inputs, double *D_eta,
		     // Monitored terms.
		     int *M_nterms, char **M_funnames, char **M_sonames,  double *M_inputs,
		     // Degree bounds.
		     int *attribs, int *maxout, int *maxin, int *minout,
		     int *minin, int *condAllDegExact, int *attriblength, 
		     // MCMC settings.
		     double *nsteps,  int *MH_interval,
		     double *burnin, double *interval,  
		     // Space for output.
		     int *collect_what,
		     double *F_sample, double *D_sample, double *M_sample,
		     int *maxedges,
		     int *newnetworktails, int *newnetworkheads, 
		     int *maxchanges,
		     int *diffnetworktime, int *diffnetworktail, int *diffnetworkhead, 
		     // Verbosity.
		     int *fVerbose,
		     int *status);

MCMCDynStatus MCMCSampleDyn(// Observed and discordant network.
			    Network *nwp,
			    // Formation terms and proposals.
			    Model *F_m, MHproposal *F_MH,
			    double *F_eta,
			    // Dissolution terms and proposals.
			    Model *D_m, MHproposal *D_MH,
			    double *D_eta,
			    // Monitored terms.
			    Model *M_m,
			    // Space for output.
			    double *F_stats, double *D_stats, double *M_stats,
			    Edge maxedges,
			    Edge maxchanges,
			    Vertex *difftime, Vertex *difftail, Vertex *diffhead,		    
			    // MCMC settings.
			    unsigned int nsteps, unsigned int MH_interval,
			    unsigned int burnin, unsigned int interval, 
			    // Verbosity.
			    int fVerbose);

MCMCDynStatus MCMCDyn1Step(// Observed and discordant network.
			   Network *nwp,
			   // Formation terms and proposals.
			   Model *F_m, MHproposal *F_MH, double *F_eta,
			   // Dissolution terms and proposals.
			   Model *D_m, MHproposal *D_MH, double *D_eta,
			   // Monitored statistics.
			   Model *M_m,
			   // Space for output.
			   unsigned log_toggles,
			   double *F_stats, double *D_stats, double *M_stats,
			   unsigned int maxchanges, Edge *nextdiffedge,
			   Vertex *difftime, Vertex *difftail, Vertex *diffhead,
			   // MCMC settings.
			   unsigned int MH_interval,
			   // Verbosity.
			   int fVerbose);
#endif
