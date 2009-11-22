#include "DynRM.h"

void MCMCDynPhase12(// Observed network.
		    int *heads, int *tails, int *n_edges,
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
		    int *fVerbose){

  Network nw[2];
  DegreeBound *bd;
  Model *F_m, *D_m;
  MHproposal F_MH, D_MH;
  DynamOrder order;
  
  
  Vertex *difftime, *diffhead, *difftail;
  difftime = (Vertex *) calloc(*maxedges,sizeof(Vertex));
  diffhead = (Vertex *) calloc(*maxedges,sizeof(Vertex));
  difftail = (Vertex *) calloc(*maxedges,sizeof(Vertex));
  
  MCMCDyn_init_common(heads, tails, *n_edges, *maxpossibleedges,
		      *n_nodes, *dflag, *bipartite, nw,
		      *order_code, &order,
		      *F_nterms, *F_funnames, *F_sonames, F_inputs, &F_m,
		      *D_nterms, *D_funnames, *D_sonames, F_inputs, &D_m,
		      attribs, maxout, maxin, minout,
		      minin, *condAllDegExact, *attriblength, &bd,
		      *F_MHproposaltype, *F_MHproposalpackage, &F_MH,
		      *D_MHproposaltype, *D_MHproposalpackage, &D_MH,
		      *fVerbose);

  MCMCSampleDynPhase12(nw, order,

		       F_m, &F_MH, theta0, 
		       init_dev, *gain,
		       *phase1n_base, *phase2n_base, *phase2sub,
		       
		       D_m, &D_MH, gamma0,
		       
		       bd,
		       *maxedges,
		       difftime, diffhead, difftail,
		       *burnin, *interval, *dyninterval,
		       *fVerbose);

  MCMCDyn_finish_common(nw, F_m, D_m, bd, &F_MH, &D_MH);

}


/*********************
 void MCMCSampleDynPhase12
*********************/
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
			  // Degree bounds.
			  DegreeBound *bd,
			  // Space for output.
			  Edge nmax,
			  Vertex *difftime, Vertex *diffhead, Vertex *difftail,
			  // MCMC settings.
			  unsigned int burnin, unsigned int interval, unsigned int dyninterval,
			  // Verbosity.
			  int fVerbose){
  int i, j;
  Edge nextdiffedge=1;
  unsigned int phase1n=phase1n_base+3*F_m->n_stats;
  double *meandev=(double*)calloc(F_m->n_stats,sizeof(double));
  
  double *ubar, *u2bar, *aDdiaginv, *D_stats;
  ubar = (double *)malloc( F_m->n_stats * sizeof(double));
  u2bar = (double *)malloc( F_m->n_stats * sizeof(double));
  aDdiaginv = (double *)malloc( F_m->n_stats * sizeof(double));
  D_stats = (double *)calloc( D_m->n_stats, sizeof(double));


  /*********************
   Burn in step. 
   *********************/
  
  if(fVerbose) Rprintf("Starting burnin of %d steps\n", burnin);
  for(i=0;i<burnin;i++)
    MCMCDyn1Step(nwp, order,
		 F_m, F_MH, theta,
		 D_m, D_MH, gamma,
		 bd,
		 0,
		 dev, D_stats,
		 nmax, &nextdiffedge,
		 difftime, diffhead, difftail,
		 dyninterval,
		 fVerbose);

  unsigned int redo, redos=3;
  do{

    redo=FALSE;

  /********************
   Phase 1
   ********************/
  
  Rprintf("Phase 1: %d steps (interval = %d)\n", phase1n,interval);
  
  for (j=0; j < F_m->n_stats; j++){
    ubar[j] = 0.0;
    u2bar[j] = 0.0;
    Rprintf("j %d %f\n",j,theta[j]);
  }

  for (i=0; i < phase1n*interval; i++){
    MCMCDyn1Step(nwp, order,
		 F_m, F_MH, theta,
		 D_m, D_MH, gamma,
		 bd,
		 0,
		 dev, D_stats,
		 nmax, &nextdiffedge,
		 difftime, diffhead, difftail,
		 dyninterval,
		 fVerbose);
    for (j=0; j<F_m->n_stats; j++){
      ubar[j]  += dev[j];
      u2bar[j] += dev[j]*dev[j];
    }
  }
  
  if (fVerbose){
    Rprintf("Returned from Phase 1\n");
    Rprintf("gain times inverse variances:\n");
  }
  
  for (j=0; j<F_m->n_stats; j++){
    aDdiaginv[j] = (u2bar[j]-ubar[j]*ubar[j]/(1.0*phase1n*interval))/(phase1n*interval);
    if( aDdiaginv[j] > 0.0){
      aDdiaginv[j] = gain/sqrt(aDdiaginv[j]);
    }else{
      aDdiaginv[j]=0.00001;
      redos--;
      if(redos>0) redo=TRUE;
    }
    if(fVerbose) Rprintf(" %f", aDdiaginv[j]);
  }
  if(fVerbose) Rprintf("\n");
  
  /********************
   Phase 2
   ********************/
  unsigned int phase2n=F_m->n_stats+7+phase2n_base;
  for(unsigned int subphase=0; subphase<phase2sub; subphase++){

    for (i=0; i < phase2n; i++){
      for(j=0; j<F_m->n_stats; j++) meandev[j]=0;
      for(j=0;j < interval;j++){
	MCMCDyn1Step(nwp, order,
		     F_m, F_MH, theta,
		     D_m, D_MH, gamma,
		     bd,
		     0,
		     dev, D_stats,
		     nmax, &nextdiffedge,
		     difftime, diffhead, difftail,
		     dyninterval,
		     fVerbose);
	for(unsigned int k=0;k<F_m->n_stats; k++)
	  meandev[k]+=dev[k];
      }
      
      /* Update theta0 */
      for (j=0; j<F_m->n_stats; j++){
	meandev[j]/=interval;
        theta[j] -= aDdiaginv[j] * meandev[j];
      }
    }
    phase2n=trunc(2.52*(phase2n-phase2n_base)+phase2n_base);
    for (j=0; j<F_m->n_stats; j++){
      aDdiaginv[j] /= 2.0;
      if (fVerbose)Rprintf("j %d theta %f ns %f\n",
			   j, theta[j], meandev[j]);
    }
    Rprintf("\n");
    
  }

  }while(redo);
  
  free(meandev);
  free(ubar);
  free(u2bar);
}
