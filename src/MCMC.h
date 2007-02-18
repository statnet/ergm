#ifndef MCMC_H
#define MCMC_H

#include "edgeTree.h"
#include "basechangeStats.h"
#include "model.h"

/*  Notes on MHproposal type:
   An MH proposal function must take two arguments:  a pointer to an 
   MHproposal structure, which holds all the information regarding the
   MH proposal; and a pointer to an array of Network structures, which 
   contain the network(s).  
   
   Each MH proposal function should check to see whether ntoggles==0
   upon being called.  If so, the MH proposal function should change
   the value of ntoggles to be the largest possible number of toggles
   required, so that this amount of memory can be allocated.
*/
typedef struct MHproposalstruct {
  void (*func)(struct MHproposalstruct*, DegreeBound*, Network*);
  Edge ntoggles;
  Vertex *togglehead;
  Vertex *toggletail;
  double ratio;
  int status;
  double *inputs; /* may be used if needed, ignored if not. */
  /* int multiplicity; Is this needed? I removed all references to
       'multiplicity' everywhere */
} MHproposal;

void MCMC_wrapper (double *heads, double *tails, double *dnedges,
		   double *dn, int *dflag, double *bipartite, 
		   int *nterms, char **funnames,
		   char **sonames, 
		   char **MHproposaltype, char **MHproposalpackage,
		   double *inputs, 
		   double *theta0, double *samplesize, double *sample,
		   double *burnin, double *interval,
		   int *newnetwork, 
		   int *fVerbose, 
		   int *attribs, int *maxout, int *maxin, int *minout,
		   int *minin, int *condAllDegExact, int *attriblength, 
		   double *maxedges,
		   double *mheads, double *mtails, double *mdnedges,
		   int *mdflag);
void MCMCSample (char *MHproposaltype, char *MHproposalpackage,
		 double *theta, double *networkstatistics, 
		 long int samplesize, long int burnin, 
		 long int interval, int hammingterm, int fVerbose,
	       	 Network *nwp, Model *m, DegreeBound *bd);
void MetropolisHastings (MHproposal *MHp,
			 double *theta, double *statistics, 
			 long int nsteps, long int *staken,
			 int hammingterm, int fVerbose,
			 Network *nwp, Model *m, DegreeBound *bd);
int CheckTogglesValid(MHproposal *MHp, DegreeBound *bd, Network *nwp);
int CheckConstrainedTogglesValid(MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MCMC_global (double *heads, double *tails, double *dnedges,
		  double *dn, int *dflag, double *bipartite,
		  int *nterms, char **funnames,
		  char **sonames, double *inputs,  double *stats); 

void MCMCDyn_wrapper (double *heads, double *tails, double *dnedges,
                   double *dn, int *dflag, double *bipartite, 
                   int *nterms, char **funnames, char **sonames, 
                   char **MHproposaltype, char **MHproposalpackage,
                   double *inputs, double *theta0, 
                   int *ndynterms, char **dynfunnames, char **dynsonames, 
                   double *dyninputs, 
		   double *samplesize, 
                   double *sample, double *burnin, double *interval,  
                   int *newnetworkhead, int *newnetworktail, 
		   int *numdissolved, int *dissolvedhead, int *dissolvedtail,
                   int *fVerbose, 
                   double *gamma, int *dyninterval,
                   int *attribs, int *maxout, int *maxin, int *minout,
                   int *minin, int *condAllDegExact, int *attriblength, 
                   double *maxedges,
                   double *mheads, double *mtails, double *mdnedges,
                   int *mdflag);
void MCMCSampleDyn (char *MHproposaltype, char *MHproposalpackage,
		 double *theta, double *networkstatistics, 
		 long int samplesize, long int burnin, 
		 long int interval, int hammingterm, int fVerbose,
	       	 double *gamma, int dyninterval,
		 Edge *numdissolve, Vertex *dissolvehead, Vertex *dissolvetail,
		 Network *nwp, Model *m, Model *mdyn, DegreeBound *bd);
void MetropolisHastingsDyn (MHproposal *MHp,
			 double *theta, double *statistics, 
			 long int nsteps, long int *staken,
			 int hammingterm, int fVerbose,
			 double *gamma, int dyninterval,
			 Edge *numdissolve,
                         Vertex *dissolvehead, Vertex *dissolvetail,
			 Network *nwp, Model *m, Model *mdyn, DegreeBound *bd);

#endif 
