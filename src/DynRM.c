#include "DynRM.h"

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
		    int *fVerbose){

  Network nw[2];
  Model *F_m, *D_m;
  MHproposal F_MH, D_MH;
  
  Vertex *difftime, *difftail, *diffhead;
  difftime = (Vertex *) calloc(*maxedges,sizeof(Vertex));
  difftail = (Vertex *) calloc(*maxedges,sizeof(Vertex));
  diffhead = (Vertex *) calloc(*maxedges,sizeof(Vertex));

  memset(newnetworktail,0,*maxedges*sizeof(int));
  memset(newnetworkhead,0,*maxedges*sizeof(int));

  MCMCDyn_init_common(tails, heads, *n_edges,
		      *n_nodes, *dflag, *bipartite, nw,
		      *F_nterms, *F_funnames, *F_sonames, F_inputs, &F_m,
		      *D_nterms, *D_funnames, *D_sonames, D_inputs, &D_m,
		      attribs, maxout, maxin, minout,
		      minin, *condAllDegExact, *attriblength,
		      *F_MHproposaltype, *F_MHproposalpackage, &F_MH,
		      *D_MHproposaltype, *D_MHproposalpackage, &D_MH,
		      *fVerbose);

  MCMCDynRMPhase2(nw,

		       F_m, F_offset, &F_MH, F_theta0, 
		       init_dev, 
		       *phase2n,
		       
		       D_m, &D_MH, D_theta0,
		       
		       *maxedges,
		       difftime, difftail, diffhead,
		       obj_history,
		       *RM_burnin, *RM_interval, *MH_interval, invGradient,
		       *fVerbose);

  /* record the final network to pass back to R */

  newnetworktail[0]=newnetworkhead[0]=EdgeTree2EdgeList(newnetworktail+1,newnetworkhead+1,nw,*maxedges);


  MCMCDyn_finish_common(nw, F_m, D_m, &F_MH, &D_MH);

}


/*********************
 void MCMCSampleDynPhase12
*********************/
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
			  int fVerbose){
  Edge nextdiffedge=1;
  double *meandev=(double*)calloc(F_m->n_stats,sizeof(double)), *meandev2=(double *)calloc( F_m->n_stats,sizeof(double)), *D_stats;
  unsigned int hist_pos=0;
  
  D_stats = (double *)calloc( D_m->n_stats, sizeof(double));

  for (unsigned int i=0; i < phase2n; i++){
    for(unsigned int j=0; j<F_m->n_stats; j++){
      meandev[j]=0;
      meandev2[j]=0;
    }
    for(unsigned int j=0;j < RM_interval;j++){
      MCMCDyn1Step(nwp,
		   F_m, F_MH, F_theta,
		   D_m, D_MH, D_theta,
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
	for(unsigned int k=0; k<F_m->n_stats; k++){
	  Rprintf("j %d F_theta %f ns %f\n",
		  k, F_theta[k], dev[k]);
            }
	Rprintf("\n");
      }
    }
    
    if(fVerbose>1){
      for (unsigned int j=0; j<F_m->n_stats; j++){
	Rprintf("j %d F_theta %f ns %f sd %f z %f\n",
            j, F_theta[j], meandev[j], sqrt(meandev2[j]-meandev[j]*meandev[j]), meandev[j]/sqrt((meandev2[j]-meandev[j]*meandev[j]))*RM_interval);
      }
      Rprintf("\n");
    }
    
    /* Record configurations and estimating equation values. */
    if(obj_history){
      for (unsigned int j=0; j<F_m->n_stats; j++){
	obj_history[hist_pos*2*F_m->n_stats+j] = F_theta[j];
	obj_history[hist_pos*2*F_m->n_stats+j+F_m->n_stats] = meandev[j]/RM_interval;
      }
      hist_pos++;
    }
    
    /* Update F_theta */
    for (unsigned int j=0; j<F_m->n_stats; j++){
      meandev[j]/=RM_interval;
      meandev2[j]/=RM_interval;
    }
    for (unsigned int j=0; j<F_m->n_stats; j++){
      for(unsigned int k=0; k<F_m->n_stats; k++)
	F_theta[j] -= invGradient[j*F_m->n_stats+k] * meandev[k]; // This may need to have k and j in invGradient switched.
    }
  }
  
  free(meandev);
  free(meandev2);
}

