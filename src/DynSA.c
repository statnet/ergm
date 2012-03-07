/*
 *  File ergm/src/DynSA.c
 *  Part of the statnet package, http://statnet.org
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) in
 *    http://statnet.org/attribution
 *
 *  Copyright 2012 the statnet development team
 */
#include "DynSA.h"

void MCMCDynSArun_wrapper(// Observed network.
			     int *tails, int *heads, int *time, int *lasttoggle, int *n_edges,
			     int *n_nodes, int *dflag, int *bipartite, 
			     // Formation terms and proposals.
			     int *F_nterms, char **F_funnames, char **F_sonames,
			     char **F_MHproposaltype, char **F_MHproposalpackage,
			     double *F_inputs, 
			     // Dissolution terms and proposals.
			     int *D_nterms, char **D_funnames, char **D_sonames, 
			     char **D_MHproposaltype, char **D_MHproposalpackage,
			     double *D_inputs,
			     // Parameter fittig.
			     double *eta0, 
			     int *M_nterms, char **M_funnames, char **M_sonames, double *M_inputs,
			     double *init_dev,
			     int *runlength,
			     double *WinvGradient,
			     double *jitter, double *dejitter,
			     // Degree bounds.
			     int *attribs, int *maxout, int *maxin, int *minout,
			     int *minin, int *condAllDegExact, int *attriblength,
			     // MCMC settings.
			     int *SA_burnin, int *SA_interval, int *MH_interval,
			     // Space for output.
			     int *maxedges, int *maxchanges,
			     int *newnetworktail, int *newnetworkhead, 
			     double *opt_history,
			     // Verbosity.
			     int *fVerbose,
			     int *status){

  Network nw[2];
  Model *F_m, *D_m, *M_m;
  MHproposal F_MH, D_MH;
  
  Vertex *difftime, *difftail, *diffhead;
  difftime = (Vertex *) calloc(*maxchanges,sizeof(Vertex));
  difftail = (Vertex *) calloc(*maxchanges,sizeof(Vertex));
  diffhead = (Vertex *) calloc(*maxchanges,sizeof(Vertex));

  memset(newnetworktail,0,*maxedges*sizeof(int));
  memset(newnetworkhead,0,*maxedges*sizeof(int));

  MCMCDyn_init_common(tails, heads, *time, lasttoggle, *n_edges,
		      *n_nodes, *dflag, *bipartite, nw,
		      *F_nterms, *F_funnames, *F_sonames, F_inputs, &F_m,
		      *D_nterms, *D_funnames, *D_sonames, D_inputs, &D_m,
		      *M_nterms, *M_funnames, *M_sonames, M_inputs, &M_m,
		      attribs, maxout, maxin, minout,
		      minin, *condAllDegExact, *attriblength,
		      *F_MHproposaltype, *F_MHproposalpackage, &F_MH,
		      *D_MHproposaltype, *D_MHproposalpackage, &D_MH,
		      *fVerbose);

  *status = MCMCDynSArun(nw,
			    
			    F_m, &F_MH,
			    D_m, &D_MH,

			    eta0, M_m,
			    init_dev, 
			    *runlength,
			    WinvGradient, jitter, dejitter,

		  *maxedges, *maxchanges,
		  difftime, difftail, diffhead,
		  opt_history,

		  *SA_burnin, *SA_interval, *MH_interval, 
		  *fVerbose);

  /* record the final network to pass back to R */

  if(*status==MCMCDyn_OK){
    newnetworktail[0]=newnetworkhead[0]=EdgeTree2EdgeList(newnetworktail+1,newnetworkhead+1,nw,*maxedges);
    *time = nw->duration_info.MCMCtimer;
    memcpy(lasttoggle, nw->duration_info.lasttoggle, sizeof(int)*(*dflag? *n_nodes*(*n_nodes-1) : (*n_nodes*(*n_nodes-1))/2));
  }

  MCMCDyn_finish_common(nw, F_m, D_m, M_m, &F_MH, &D_MH);
  free(difftime);
  free(difftail);
  free(diffhead);
}


/*********************
 void MCMCSampleDynPhase12
*********************/
MCMCDynStatus MCMCDynSArun(// Observed and discordant network.
			      Network *nwp,
			      // Formation terms and proposals.
			      Model *F_m, MHproposal *F_MH,
			      // Dissolution terms and proposals.
			      Model *D_m, MHproposal *D_MH,
			      // Model fitting.
			      double *eta, 
			      Model *M_m,
			      double *dev, // DEViation of the current network's targeted statistics from the target statistics.
			      int runlength,
			      double *WinvGradient, double *jitter, double *dejitter,
			      
			      // Space for output.
			      Edge maxedges, Edge maxchanges,
			      Vertex *difftime, Vertex *difftail, Vertex *diffhead,
			      double *opt_history,
			      // MCMC settings.
			      unsigned int SA_burnin, unsigned int SA_interval, unsigned int MH_interval,
			      // Verbosity.
			      int fVerbose){
  Edge nextdiffedge=1;
  unsigned int hist_pos=0, p=F_m->n_stats+D_m->n_stats, n, rowsize = p*2 + M_m->n_stats;
  double *meandev=(double*)calloc(M_m->n_stats,sizeof(double)), *last_jitter=(double*)calloc(p,sizeof(double));


  for (unsigned int i=0; i < runlength; i++){
    n = 0;
    for(unsigned int j=0; j<M_m->n_stats; j++){
      meandev[j]=0;
    }

    // Jitter parameters
    for (unsigned int j=0; j<p; j++){
      if(jitter[j]!=0){
	last_jitter[j] = rnorm(0,jitter[j]);
	eta[j] += last_jitter[j];
      }
    }

    for(unsigned int j=0;j < SA_interval;j++){
      MCMCDynStatus status = MCMCDyn1Step(nwp,
					  F_m, F_MH, eta,
					  D_m, D_MH, eta+F_m->n_stats,
					  M_m,
					  0,
					  NULL, NULL, dev,
					  maxchanges, &nextdiffedge,
					  difftime, difftail, diffhead,
					  MH_interval,
					  fVerbose);
      if(status==MCMCDyn_TOO_MANY_CHANGES)
        return MCMCDyn_TOO_MANY_CHANGES;
      
      if(nwp->nedges >= maxedges-1)
	return MCMCDyn_TOO_MANY_EDGES;

      for(unsigned int k=0;k<M_m->n_stats; k++){
	meandev[k]+=dev[k]*j;
	n+=j;
      }
      if (fVerbose>2){
	for(unsigned int k=0; k<p; k++){
	  Rprintf("eta[%d] = %f\n", k, eta[k]);
	}
	for(unsigned int k=0; k<M_m->n_stats; k++){
	  Rprintf("M_dev[%d] = %f\n", k, dev[k]);
	}

	Rprintf("\n");
      }
    }
    
    if(fVerbose>1){
      for(unsigned int k=0; k<p; k++){
	Rprintf("eta[%d] = %f\n", k, eta[k]);
      }
      for(unsigned int k=0; k<M_m->n_stats; k++){
	Rprintf("meandev[%d] = %f\n", k, meandev[k]/n);
      }
      
      Rprintf("\n");
    }

    // Evaluate mean deviations.
    for (unsigned int j=0; j<M_m->n_stats; j++){
      meandev[j]/=n;
    }
    
    // Record configurations and estimating equation values.
 
    
    for(unsigned int j=0; j<p; j++){
      opt_history[hist_pos*rowsize+j] = eta[j];
    }
    for(unsigned int j=0; j<p; j++){
      opt_history[hist_pos*rowsize+p+j] = last_jitter[j];
    }
    for(unsigned int j=0; j<M_m->n_stats; j++){
      opt_history[hist_pos*rowsize+p+p+j] = meandev[j];
    }
    hist_pos++;
  
    // Update formation and dissolution parameters, and cancel the effect of jitter.
    // eta[t+1] = eta[t] - a*(G^-1)*W*(d[t] - G*jit[t])
    //          = eta[t] - a*(G^-1)*W*d[t] + a*(G^-1)*W*G*jit[t]
    for(unsigned int k=0; k<M_m->n_stats; k++){
      for (unsigned int j=0; j<p; j++){
	eta[j] -= WinvGradient[k*p+j] * meandev[k];
      }
    }
    for(unsigned int k=0; k<p; k++){
      for (unsigned int j=0; j<p; j++){
	eta[j] += dejitter[k*p+j] * last_jitter[k];
      }
    }

    // Undo jitter
    for (unsigned int j=0; j<p; j++){
      if(jitter[j]!=0){
	eta[j] -= last_jitter[j];
      }
    }
  }

  free(meandev);
  free(last_jitter);
  return MCMCDyn_OK;
}

