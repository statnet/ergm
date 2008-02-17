#ifndef MCMC_H
#define MCMC_H

#include "edgetree.h"
#include "changestats.h"
#include "model.h"

/* Macros to test for logical inequality (XOR) and logical equality (XNOR). */
#define XOR(a,b) (((a)==0) != ((b)==0))
#define XNOR(a,b) (((a)==0) == ((b)==0))

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

void MCMC_wrapper (int *heads, int *tails, int *dnedges,
                   int *maxpossibleedges,
                   int *dn, int *dflag, int *bipartite, 
                   int *nterms, char **funnames,
                   char **sonames, 
                   char **MHproposaltype, char **MHproposalpackage,
                   double *inputs, double *theta0, int *samplesize, 
                   double *sample, int *burnin, int *interval,  
                   int *newnetworkheads, 
                   int *newnetworktails, 
                   int *fVerbose, 
                   int *attribs, int *maxout, int *maxin, int *minout,
                   int *minin, int *condAllDegExact, int *attriblength, 
                   int *maxedges,
                   int *mheads, int *mtails, int *mdnedges);
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
void MCMC_global (int *heads, int *tails, int *dnedges,
      int *maxpossibleedges,
		  int *dn, int *dflag,  int *bipartite,
		  int *nterms, char **funnames,
		  char **sonames, double *inputs,  double *stats);
void MCMCPhase12 (int *heads, int *tails, int *dnedges,
      int *maxpossibleedges,
		  int *dn, int *dflag, int *bipartite, 
		  int *nterms, char **funnames,
		  char **sonames, 
		  char **MHproposaltype, char **MHproposalpackage,
		  double *inputs, 
		  double *theta0, int *samplesize,
		  double *gain, double *meanstats, int *phase1, int *nsub,
		  double *sample, int *burnin, int *interval,  
		  int *newnetworkheads, 
		  int *newnetworktails, 
		  int *fVerbose, 
		  int *attribs, int *maxout, int *maxin, int *minout,
		  int *minin, int *condAllDegExact, int *attriblength, 
		  int *maxedges,
		  int *mheads, int *mtails, int *mdnedges);

void MCMCSamplePhase12 (char *MHproposaltype, char *MHproposalpackage,
  double *theta, double gain, double *meanstats,
  int nphase1, int nsubphases, double *networkstatistics, 
  long int samplesize, long int burnin, 
  long int interval, int hammingterm, int fVerbose,
  Network *nwp, Model *m, DegreeBound *bd);

void MH_init(MHproposal *MH, 
	     char *MHproposaltype, char *MHproposalpackage, 
	     int fVerbose,
	     Network *nwp, DegreeBound *bd);

void MH_free(MHproposal *MH);

void ChangeStats(unsigned int ntoggles, Vertex *togglehead, Vertex *toggletail, Network *nwp, Model *m);

#endif
