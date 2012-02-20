#include "DynRM.h"

void MCMCDynRMPhase2_wrapper(// Observed network.
		    int *tails, int *heads, int *n_edges,
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
		    double *invGradient,
		    double *jitter,
		    // Degree bounds.
		    int *attribs, int *maxout, int *maxin, int *minout,
		    int *minin, int *condAllDegExact, int *attriblength,
		    // MCMC settings.
		    int *RM_burnin, int *RM_interval, int *MH_interval,
		    // Space for output.
		    int *maxedges,
		    int *newnetworktail, int *newnetworkhead, 
		    double *opt_history,
		    // Verbosity.
		    int *fVerbose){

  Network nw[2];
  Model *F_m, *D_m, *M_m;
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
		      *M_nterms, *M_funnames, *M_sonames, M_inputs, &M_m,
		      attribs, maxout, maxin, minout,
		      minin, *condAllDegExact, *attriblength,
		      *F_MHproposaltype, *F_MHproposalpackage, &F_MH,
		      *D_MHproposaltype, *D_MHproposalpackage, &D_MH,
		      *fVerbose);

  MCMCDynRMPhase2(nw,

		  F_m, &F_MH, F_eta0,  
		  D_m, &D_MH, D_eta0,

		  M_m,
		  init_dev, 
		  *phase2n,
		  invGradient, jitter,

		  *maxedges,
		  difftime, difftail, diffhead,
		  opt_history,

		  *RM_burnin, *RM_interval, *MH_interval, 
		  *fVerbose);

  /* record the final network to pass back to R */

  newnetworktail[0]=newnetworkhead[0]=EdgeTree2EdgeList(newnetworktail+1,newnetworkhead+1,nw,*maxedges);


  MCMCDyn_finish_common(nw, F_m, D_m, M_m, &F_MH, &D_MH);

}


/*********************
 void MCMCSampleDynPhase12
*********************/
void MCMCDynRMPhase2(// Observed and discordant network.
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
			  double *invGradient, double *jitter,

			  // Space for output.
			  Edge nmax,
			  Vertex *difftime, Vertex *difftail, Vertex *diffhead,
			  double *opt_history,
			  // MCMC settings.
			  unsigned int RM_burnin, unsigned int RM_interval, unsigned int MH_interval,
			  // Verbosity.
			  int fVerbose){
  Edge nextdiffedge=1;
  double *meandev=(double*)calloc(M_m->n_stats,sizeof(double)), *meandev2=(double *)calloc( M_m->n_stats,sizeof(double)), *last_jitter=(double*)calloc(F_m->n_stats+D_m->n_stats,sizeof(double));
  unsigned int hist_pos=0;
  

  for (unsigned int i=0; i < phase2n; i++){
    for(unsigned int j=0; j<M_m->n_stats; j++){
      meandev[j]=0;
      meandev2[j]=0;
    }
    for(unsigned int j=0;j < RM_interval;j++){
      MCMCDyn1Step(nwp,
		   F_m, F_MH, F_eta,
		   D_m, D_MH, D_eta,
		   M_m,
		   0,
		   NULL, NULL, dev,
		   nmax, &nextdiffedge,
		   difftime, difftail, diffhead,
		   MH_interval,
		   fVerbose);
      for(unsigned int k=0;k<M_m->n_stats; k++){
	meandev[k]+=dev[k];
	meandev2[k]+=dev[k]*dev[k];
      }
      if (fVerbose>2){
	for(unsigned int k=0; k<F_m->n_stats; k++){
	  Rprintf("F_eta[%d] = %f\n", k, F_eta[k]);
	}
	for(unsigned int k=0; k<D_m->n_stats; k++){
	  Rprintf("D_eta[%d] = %f\n", k, D_eta[k]);
	}
	for(unsigned int k=0; k<M_m->n_stats; k++){
	  Rprintf("M_dev[%d] = %f\n", k, dev[k]);
	}

	Rprintf("\n");
      }
    }
    
    if(fVerbose>1){
      for(unsigned int k=0; k<F_m->n_stats; k++){
	Rprintf("F_eta[%d] = %f\n", k, F_eta[k]);
      }
      for(unsigned int k=0; k<D_m->n_stats; k++){
	Rprintf("D_eta[%d] = %f\n", k, D_eta[k]);
      }
      for(unsigned int k=0; k<M_m->n_stats; k++){
	Rprintf("meandev[%d] = %f, sd[%d] = %f, z[%d] = %f\n", k, meandev[k], k, sqrt(meandev2[k]-meandev[k]*meandev[k]), k, meandev[k]/sqrt((meandev2[k]-meandev[k]*meandev[k]))*RM_interval);
      }
      
      Rprintf("\n");
    }

    // Evaluate mean deviations.
    for (unsigned int j=0; j<M_m->n_stats; j++){
      meandev[j]/=RM_interval;
      meandev2[j]/=RM_interval;
    }

    
    // Record configurations and estimating equation values.
    unsigned int rowsize = F_m->n_stats + D_m->n_stats + M_m->n_stats;
    
    for(unsigned int j=0; j<F_m->n_stats; j++){
      opt_history[hist_pos*rowsize+j] = F_eta[j];
    }
    for(unsigned int j=0; j<D_m->n_stats; j++){
      opt_history[hist_pos*rowsize+F_m->n_stats+j] = D_eta[j];
    }
    for(unsigned int j=0; j<M_m->n_stats; j++){
      opt_history[hist_pos*rowsize+F_m->n_stats+D_m->n_stats+j] = meandev[j];
    }
    hist_pos++;
  
    // Update formation and dissolution parameters, with optional jitter.
    for(unsigned int k=0; k<M_m->n_stats; k++){
      for (unsigned int j=0; j<F_m->n_stats; j++){
	F_eta[j] -= invGradient[j*M_m->n_stats+k] * meandev[k];
      }
      for (unsigned int j=0; j<D_m->n_stats; j++){
	D_eta[j] -= invGradient[(F_m->n_stats+j)*M_m->n_stats+k] * meandev[k];
      }
    }

    for (unsigned int j=0; j<F_m->n_stats; j++){
      if(jitter[j]!=0){
	F_eta[j] -= last_jitter[j];
	last_jitter[j] = rnorm(0,jitter[j]);
	F_eta[j] += last_jitter[j];
      }
    }
    for (unsigned int j=0; j<D_m->n_stats; j++){
      if(jitter[F_m->n_stats+j]!=0){
	D_eta[j] -= last_jitter[F_m->n_stats+j];
	last_jitter[F_m->n_stats+j] = rnorm(0,jitter[F_m->n_stats+j]);
	D_eta[j]+=rnorm(0,jitter[F_m->n_stats+j]);
      }
    }
  }
  
  free(meandev);
  free(meandev2);
  free(last_jitter);
}

