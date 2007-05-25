#ifndef MCMCDYN_H
#define MCMCDYN_H

typedef enum {
  DissThenForm=1,
  DissAndForm=2,
  FormThenDiss=3,
  // These modes are used in debugging:
  FormOnly=4,
  DissOnly=5
} DynamOrder;

void MCMCDyn_wrapper(// Starting network.
		     int *heads, int *tails, int *dnedges,
		     int *dn, int *dflag, int *bipartite,
		     // Ordering of formation and dissolution.
		     int *order_code,
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
		     double *nsteps,  int *dyninterval,
		     double *burnin, double *interval,  
		     // Space for output.
		     double *F_sample, double *D_sample, 
		     int *newnetworkhead, int *newnetworktail, 
		     double *maxedges,
		     int *diffnetworktime, int *diffnetworkhead, int *diffnetworktail, 
		     // Verbosity.
		     int *fVerbose);

void MCMCSampleDyn(// Observed and discordant network.
		   Network *nwp,
		   // Ordering of formation and dissolution.
		   DynamOrder order,
		   // Formation terms and proposals.
		   Model *F_m, MHproposal *F_MH, double *theta,
		   // Dissolution terms and proposals.
		   Model *D_m, MHproposal *D_MH, double *gamma,
		   // Degree bounds.
		   DegreeBound *bd,
		   // Space for output.
		   double *F_stats, double *D_stats,// Do we still need these?
		   Edge nmax,
		   Vertex *difftime, Vertex *diffhead, Vertex *difftail,		    
		   // MCMC settings.
		   unsigned int nsteps, unsigned int dyninterval,
		   unsigned int burnin, unsigned int interval, 
		   // Verbosity.
		   int fVerbose);

void MCMCDyn1Step(Network *nwp,
		  DynamOrder order,
		  Model *F_m, MHproposal *F_MH, double *theta,
		  Model *D_m, MHproposal *D_MH, double *gamma,
		  DegreeBound *bd,
		  unsigned log_toggles,
		  double *F_stats, double *D_stats,
		  unsigned int nmax, unsigned int *nextdiffedge,
		  Vertex *difftime, Vertex *diffhead, Vertex *difftail,
		  unsigned int dyninterval,
		  int fVerbose);

void MCMCDynPhase12(// Observed network.
		    int *heads, int *tails, int *n_edges,
		    int *n_nodes, int *dflag, int *bipartite, 
		    // Ordering of formation and dissolution.
		    int *order_code, 
		    // Formation terms and proposals.
		    int *F_nterms, char **F_funnames, char **F_sonames, 
		    char **F_MHproposaltype, char **F_MHproposalpackage,
		    double *F_inputs, double *theta0, 
		    // Formation parameter fitting.
		    double *gain, double *meanstats, 
		    int *phase1, int *nsub,
		    // Dissolution terms and proposals.
		    int *D_nterms, char **D_funnames, char **D_sonames, 
		    char **D_MHproposaltype, char **D_MHproposalpackage,
		    double *D_inputs, double *gamma0,
		    // Dissolution parameter fitting --- to add later? -PK
		    // Degree bounds.
		    int *attribs, int *maxout, int *maxin, int *minout,
		    int *minin, int *condAllDegExact, int *attriblength, 
		    // MCMC settings.
		    double *nsteps, int *dyninterval,
		    double *burnin, double *interval, 
		    // Space for output.
		    double *F_sample, double *D_sample, 
		    int *newnetworkhead, int *newnetworktail, 
		    double *maxedges,
		    int *diffnetworktime, int *diffnetworkhead, int *diffnetworktail, 
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
			  double gain, double *meanstats, int nphase1, int nsubphases,
			  // Dissolution terms and proposals.
			  Model *D_m, MHproposal *D_MH,
			  double *gamma, 
			  // Dissolution parameter fitting --- to add later? -PK
			  // Degree bounds.
			  DegreeBound *bd,
			  // Space for output.
			  double *F_stats,  double *D_stats, // Do we still need theese?
			  Edge nmax,
			  Vertex *difftime, Vertex *diffhead, Vertex *difftail,
			  // MCMC settings.
			  unsigned int samplesize, unsigned int burnin, 
			  unsigned int interval, unsigned int dyninterval,
			  // Verbosity.
			  int fVerbose);

#endif
