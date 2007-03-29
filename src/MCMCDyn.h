#ifndef MCMCDYN_H
#define MCMCDYN_H

typedef enum {DissThenForm=1, DissAndForm=2, FormThenDiss=3} DynamOrder;


void MCMCDyn_wrapper (int *order_code, double *heads, double *tails, double *dnedges,
                   double *dn, int *dflag, double *bipartite, 
                   int *nterms, char **funnames, char **sonames, 
                   char **MHproposaltype, char **MHproposalpackage,
                   double *inputs, double *theta0, 
                   int *ndynterms, char **dynfunnames, char **dynsonames, 
                   double *dyninputs, 
		   double *samplesize, 
                   double *sample, double *burnin, double *interval,  
                   int *newnetworkhead, int *newnetworktail, 
                   int *diffnetworktime, int *diffnetworkhead, int *diffnetworktail, 
                   int *dissnetworktime, int *dissnetworkhead, int *dissnetworktail, 
                   int *fVerbose, 
                   double *gamma, int *dyninterval,
                   int *attribs, int *maxout, int *maxin, int *minout,
                   int *minin, int *condAllDegExact, int *attriblength, 
                   double *maxedges,
                   double *mheads, double *mtails, double *mdnedges,
                   int *mdflag);
void MCMCSampleDyn (DynamOrder order, char *MHproposaltype, char *MHproposalpackage,
		 double *theta, double *networkstatistics, 
		 long int samplesize, long int burnin, 
		 long int interval, int hammingterm, int fVerbose,
	       	 double *gamma, int dyninterval,
		 Edge *nmax,
		 Vertex *dissolvetime, Vertex *dissolvehead, Vertex *dissolvetail,
		 Vertex *disstime, Vertex *disshead, Vertex *disstail,
		 Vertex *difftime, Vertex *diffhead, Vertex *difftail,
		 Network *nwp, Model *m, Model *mdyn, DegreeBound *bd);
void MetropolisHastingsDyn (DynamOrder order, MHproposal *MHp,
		 double *theta, double *statistics, 
		 long int nsteps, long int *staken,
		 int hammingterm, int fVerbose,
		 double *gamma, int dyninterval,
		 Edge *nmax,
                 Vertex *dissolvetime, Vertex *dissolvehead, Vertex *dissolvetail,
		 Vertex *disstime, Vertex *disshead, Vertex *disstail,
		 Vertex *difftime, Vertex *diffhead, Vertex *difftail,
		 Network *nwp, Model *m, Model *mdyn, DegreeBound *bd);

void MCMCDynPhase2 (int *order_code, double *heads, double *tails, double *dnedges,
                   double *dn, int *dflag, double *bipartite, 
                   int *nterms, char **funnames, char **sonames, 
                   char **MHproposaltype, char **MHproposalpackage,
                   double *inputs, 
		   double *theta0, double *aDdiaginv, 
                   int *ndynterms, char **dynfunnames, char **dynsonames, 
                   double *dyninputs, 
		   double *samplesize, 
                   double *sample, double *burnin, double *interval,  
                   int *newnetworkhead, int *newnetworktail, 
                   int *diffnetworktime, int *diffnetworkhead, int *diffnetworktail, 
                   int *dissnetworktime, int *dissnetworkhead, int *dissnetworktail, 
                   int *fVerbose, 
                   double *gamma, int *dyninterval,
                   int *attribs, int *maxout, int *maxin, int *minout,
                   int *minin, int *condAllDegExact, int *attriblength, 
                   double *maxedges,
                   double *mheads, double *mtails, double *mdnedges,
                   int *mdflag);

void MCMCSampleDynPhase2 (DynamOrder order, char *MHproposaltype, char *MHproposalpackage,
  double *theta, double *aDdiaginv, double *networkstatistics, 
  long int samplesize, long int burnin, 
  long int interval, int hammingterm, int fVerbose,
  double *gamma, int dyninterval,
  Edge *nmax,
  Vertex *dissolvetime, Vertex *dissolvehead, Vertex *dissolvetail,
  Vertex *disstime, Vertex *disshead, Vertex *disstail,
  Vertex *difftime, Vertex *diffhead, Vertex *difftail,
  Network *nwp, Model *m, Model *mdyn, DegreeBound *bd);

#endif
