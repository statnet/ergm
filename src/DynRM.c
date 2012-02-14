#include "DynRM.h"

void MCMCDynPhase12(// Observed network.
		    int *tails, int *heads, int *n_edges,
		    int *maxpossibleedges,
		    int *n_nodes, int *dflag, int *bipartite, 
		    // Formation terms and proposals.
		    int *F_nterms, char **F_funnames, char **F_sonames, int *F_offset,
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
		    int *RM_burnin, int *RM_interval, int *MH_interval,
		    // Space for output.
		    int *maxedges,
		    // Verbosity.
		    int *fVerbose){

  Network nw[2];
  Model *F_m, *D_m;
  MHproposal F_MH, D_MH;
  
  Vertex *difftime, *difftail, *diffhead;
  difftime = (Vertex *) calloc(*maxedges,sizeof(Vertex));
  difftail = (Vertex *) calloc(*maxedges,sizeof(Vertex));
  diffhead = (Vertex *) calloc(*maxedges,sizeof(Vertex));
  
  MCMCDyn_init_common(tails, heads, *n_edges, *maxpossibleedges,
		      *n_nodes, *dflag, *bipartite, nw,
		      *F_nterms, *F_funnames, *F_sonames, F_inputs, &F_m,
		      *D_nterms, *D_funnames, *D_sonames, D_inputs, &D_m,
		      attribs, maxout, maxin, minout,
		      minin, *condAllDegExact, *attriblength,
		      *F_MHproposaltype, *F_MHproposalpackage, &F_MH,
		      *D_MHproposaltype, *D_MHproposalpackage, &D_MH,
		      *fVerbose);

  MCMCSampleDynPhase12(nw,

		       F_m, F_offset, &F_MH, theta0, 
		       init_dev, *gain,
		       *phase1n_base, *phase2n_base, *phase2sub,
		       
		       D_m, &D_MH, gamma0,
		       
		       *maxedges,
		       difftime, difftail, diffhead,
		       *RM_burnin, *RM_interval, *MH_interval,
		       *fVerbose);

  MCMCDyn_finish_common(nw, F_m, D_m, &F_MH, &D_MH);

}


/*********************
 void MCMCSampleDynPhase12
*********************/
void MCMCSampleDynPhase12(// Observed and discordant network.
			  Network *nwp,
			  // Formation terms and proposals.
			  Model *F_m, int *F_offset, MHproposal *F_MH,
			  double *theta, 
			  // Formation parameter fitting.
			  double *dev, // DEViation of the current network's formation statistics from the target statistics.
			  double gain,
			  int phase1n_base, int phase2n_base, int phase2sub,
			  // Dissolution terms and proposals.
			  Model *D_m, MHproposal *D_MH,
			  double *gamma, 
			  // Dissolution parameter fitting --- to add later? -PK
			  // Space for output.
			  Edge nmax,
			  Vertex *difftime, Vertex *difftail, Vertex *diffhead,
			  // MCMC settings.
			  unsigned int RM_burnin, unsigned int RM_interval, unsigned int MH_interval,
			  // Verbosity.
			  int fVerbose){
  int i, j;
  Edge nextdiffedge=1;
  unsigned int phase1n=phase1n_base+3*F_m->n_stats;
  double *meandev=(double*)calloc(F_m->n_stats,sizeof(double)), *meandev2=(double *)calloc( F_m->n_stats,sizeof(double)), *aDdiaginv, *D_stats;
  
  aDdiaginv = (double *)calloc( F_m->n_stats, sizeof(double));
  D_stats = (double *)calloc( D_m->n_stats, sizeof(double));


  /*********************
   Burn in step. 
   *********************/
  
  if(fVerbose) Rprintf("Starting burnin of %d steps\n", RM_burnin);
  for(i=0;i<RM_burnin;i++)
    MCMCDyn1Step(nwp,
		 F_m, F_MH, theta,
		 D_m, D_MH, gamma,
		 0,
		 dev, D_stats,
		 nmax, &nextdiffedge,
		 difftime, difftail, diffhead,
		 MH_interval,
		 fVerbose);

  unsigned int redo, redos=3;
  do{

    redo=FALSE;

    /********************
    Phase 1
    ********************/
  
    Rprintf("Phase 1: %d steps (RM_interval = %d)\n", phase1n,RM_interval);
  
    for (j=0; j < F_m->n_stats; j++){
      meandev[j] = 0.0;
      meandev2[j] = 0.0;
      Rprintf("j %d %f\n",j,theta[j]);
    }

    for (i=0; i < phase1n*RM_interval; i++){
      MCMCDyn1Step(nwp,
        F_m, F_MH, theta,
        D_m, D_MH, gamma,
        0,
        dev, D_stats,
        nmax, &nextdiffedge,
        difftime, difftail, diffhead,
        MH_interval,
        fVerbose);
      for (j=0; j<F_m->n_stats; j++){
        meandev[j]  += dev[j];
        meandev2[j] += dev[j]*dev[j];
      }
    }
  
    if (fVerbose){
      Rprintf("Returned from Phase 1\n");
      Rprintf("gain times inverse variances:\n");
    }
  
    for (j=0; j<F_m->n_stats; j++){
      aDdiaginv[j] = (meandev2[j]-meandev[j]*meandev[j]/(1.0*phase1n*RM_interval))/(phase1n*RM_interval);
      if(!F_offset[j]){
	if( aDdiaginv[j] > 0.0){
	  aDdiaginv[j] = gain/sqrt(aDdiaginv[j]);
	}else{
	  aDdiaginv[j]=0.0001;
	  redos--;
	  if(redos>0) redo=TRUE;
	}
      }
      if(fVerbose) Rprintf(" %f", aDdiaginv[j]);
    }
    if(fVerbose) Rprintf("\n");

    if(fVerbose){
      Rprintf("theta0 statistics:");
      for (j=0; j<F_m->n_stats; j++){
        Rprintf("%f ",  meandev[j]/(phase1n*RM_interval));
      }
      Rprintf("\n");
    }
    
    /********************
    Phase 2
    ********************/
    unsigned int phase2n=F_m->n_stats+7+phase2n_base;
    for(unsigned int subphase=0; subphase<phase2sub; subphase++){
      
      for (i=0; i < phase2n; i++){
        for(j=0; j<F_m->n_stats; j++){
          meandev[j]=0;
          meandev2[j]=0;
        }
        for(j=0;j < RM_interval;j++){
          MCMCDyn1Step(nwp,
            F_m, F_MH, theta,
            D_m, D_MH, gamma,
            0,
            dev, D_stats,
            nmax, &nextdiffedge,
            difftime, difftail, diffhead,
            MH_interval,
            fVerbose);
          for(unsigned int k=0;k<F_m->n_stats; k++){
            meandev[k]+=dev[k];
            meandev2[k]+=dev[k]*dev[k];
          }
          if (fVerbose>2){
            for (unsigned int k=0; k<F_m->n_stats; k++){
              Rprintf("j %d theta %f ns %f\n",
              k, theta[k], dev[k]);
            }
            Rprintf("\n");
          }
        }
        
        if(fVerbose>1){
          for (j=0; j<F_m->n_stats; j++){
            Rprintf("j %d theta %f ns %f sd %f z %f\n",
            j, theta[j], meandev[j], sqrt(meandev2[j]-meandev[j]*meandev[j]), meandev[j]/sqrt((meandev2[j]-meandev[j]*meandev[j]))*RM_interval);
          }
          Rprintf("\n");
        }
        
        /* Update theta0 */
        for (j=0; j<F_m->n_stats; j++){
          meandev[j]/=RM_interval;
          meandev2[j]/=RM_interval;
          theta[j] -= aDdiaginv[j] * meandev[j];
        }
      }
      
      phase2n=trunc(2.52*(phase2n-phase2n_base)+phase2n_base);
      for (j=0; j<F_m->n_stats; j++){
        aDdiaginv[j] /= 2.0;
        if (fVerbose)Rprintf("j %d theta %f ns %f sd %f z %f\n",
          j, theta[j], meandev[j], sqrt(meandev2[j]-meandev[j]*meandev[j]), meandev[j]/sqrt((meandev2[j]-meandev[j]*meandev[j]))*RM_interval);
      }
      Rprintf("\n");
    }
    
  }while(redo);
  
  free(meandev);
  free(meandev2);
}

