#include "DynSPSA.h"

void MCMCDynSPSA_wrapper(// Observed network.
		    int *tails, int *heads, int *n_edges,
		    int *maxpossibleedges,
		    int *n_nodes, int *dflag, int *bipartite, 
		    // Ordering of formation and dissolution.
		    int *order_code, 
		    // Formation terms and proposals.
		    int *F_nterms, char **F_funnames, char **F_sonames, 
		    char **F_MHproposaltype, char **F_MHproposalpackage,
		    double *F_inputs, double *F_theta0, 
		    double *init_dev,
		    // SPSA settings.
		    double *a,
		    double *alpha,
		    double *A,
		    double *c,
		    double *gamma,
		    int *iterations,
		    // Dissolution terms and proposals.
		    int *D_nterms, char **D_funnames, char **D_sonames, 
		    char **D_MHproposaltype, char **D_MHproposalpackage,
		    double *D_inputs, double *D_theta0,
		    // Degree bounds.
		    int *attribs, int *maxout, int *maxin, int *minout,
		    int *minin, int *condAllDegExact, int *attriblength,
		    // MCMC settings.
		    int *burnin, int *interval, int *dyninterval,
		    // Space for output.
		    int *maxedges,
		    double *obj_history,
		    // Verbosity.
		    int *fVerbose){

  Network nw[2];
  Model *F_m, *D_m;
  MHproposal F_MH, D_MH;
  DynamOrder order;
  
  
  Vertex *difftime, *difftail, *diffhead;
  difftime = (Vertex *) calloc(*maxedges,sizeof(Vertex));
  difftail = (Vertex *) calloc(*maxedges,sizeof(Vertex));
  diffhead = (Vertex *) calloc(*maxedges,sizeof(Vertex));
  
  MCMCDyn_init_common(tails, heads, *n_edges, *maxpossibleedges,
		      *n_nodes, *dflag, *bipartite, nw,
		      *order_code, &order,
		      *F_nterms, *F_funnames, *F_sonames, F_inputs, &F_m,
		      *D_nterms, *D_funnames, *D_sonames, D_inputs, &D_m,
		      attribs, maxout, maxin, minout,
		      minin, *condAllDegExact, *attriblength,
		      *F_MHproposaltype, *F_MHproposalpackage, &F_MH,
		      *D_MHproposaltype, *D_MHproposalpackage, &D_MH,
		      *fVerbose);

  MCMCDynSPSA(nw, order,
	      
	      F_m, &F_MH, F_theta0, 
	      init_dev, *a, *alpha, *A, *c, *gamma, *iterations,
	      
	      D_m, &D_MH, D_theta0,
	      
	      *maxedges,
	      obj_history,
	      difftime, difftail, diffhead,
	      *burnin, *interval, *dyninterval,
	      *fVerbose);
  
  MCMCDyn_finish_common(nw, F_m, D_m, &F_MH, &D_MH);

}


double MCMCSampleDynObjective(Network *nwp,
			      // Ordering of formation and dissolution.
			      DynamOrder order,
			      // Formation terms and proposals.
			      Model *F_m, MHproposal *F_MH, double *F_theta,
			      // Dissolution terms and proposals.
			      Model *D_m, MHproposal *D_MH, double *D_theta,
			      double *F_stats, double *D_stats,
			      unsigned int nmax,
			      Vertex *difftime, Vertex *difftail, Vertex *diffhead,
			      // MCMC settings.
			      unsigned int dyninterval,
			      unsigned int burnin,
			      unsigned int S,
			      double *F_stats_acc, double* F_stats2_acc, int *use_var,
			      int fVerbose){

  // These are the same... for now.
  unsigned int n_stats=F_m->n_stats;
  unsigned int n_par=F_m->n_stats;

  if(fVerbose){
    Rprintf("MCMC Run: %d steps after discarding %d\n F_theta=[ ",S,burnin);
    for(unsigned int k=0; k<n_par; k++) Rprintf("%f ",F_theta[k]);
    Rprintf("]\n");
  }

  memset(F_stats_acc,0,sizeof(double)*n_stats);			
  if(F_stats2_acc) memset(F_stats2_acc,0,sizeof(double)*n_stats);			
  for(unsigned int s=0; s<burnin; s++)					
    MCMCDyn1Step(nwp, order, F_m, F_MH, F_theta, D_m, D_MH, D_theta, 0, F_stats, D_stats, nmax, NULL, difftime, difftail, diffhead, dyninterval, 0); 
  for(unsigned int s=0; s<S; s++){					
    MCMCDyn1Step(nwp, order, F_m, F_MH, F_theta, D_m, D_MH, D_theta, 0, F_stats, D_stats, nmax, NULL, difftime, difftail, diffhead, dyninterval, 0); 
    for(unsigned int k=0; k<n_stats; k++){				
      F_stats_acc[k]+=F_stats[k];					
      if(F_stats2_acc){
	F_stats2_acc[k]+=F_stats[k]*F_stats[k];
      }
    }
  }															
  double result=0;
  unsigned int var_good=TRUE;
  for(unsigned int k=0; k<n_stats; k++){
    double var;
    if(F_stats2_acc){
      var=fmax(F_stats2_acc[k]-F_stats_acc[k]*F_stats_acc[k]/S,0)/S;
      
      if(var/(F_stats2_acc[k]+F_stats_acc[k]*F_stats_acc[k]/S)<0.0001) var_good=FALSE;
    }else var=1;
    if(use_var && *use_var<0) var=1;
    
    result+=F_stats_acc[k]*F_stats_acc[k]/S/S/var;
  }

  if(var_good && use_var && *use_var<1) (*use_var)++;

  if(fVerbose){
    Rprintf(" g(Y)=[ ");
    for(unsigned int k=0; k<n_stats; k++) Rprintf("%f(%f) ",F_stats_acc[k]/S,sqrt((F_stats2_acc[k]-F_stats_acc[k]*F_stats_acc[k]/S)/S));
    Rprintf("]\n Objective=%f\n",result);
  }

  return result;
}


/*********************
 void MCMCDynSPSA
*********************/
void MCMCDynSPSA(// Observed and discordant network.
		       Network *nwp,
		       // Ordering of formation and dissolution.
		       DynamOrder order,
		       // Formation terms and proposals.
		       Model *F_m, MHproposal *F_MH,
		       double *F_theta, 
		       // Formation parameter fitting.
		       double *dev, // DEViation of the current network's formation statistics from the target statistics.
		       // SPSA settings.
		       double a,
		       double alpha,
		       double A,
		       double c,
		       double gamma,
		       unsigned int iterations,
		       // Dissolution terms and proposals.
		       Model *D_m, MHproposal *D_MH,
		       double *D_theta, 
		       // Dissolution parameter fitting --- to add later? -PK
		       // Space for output.
		       Edge nmax,
		       double *obj_history,
		       Vertex *difftime, Vertex *difftail, Vertex *diffhead,
		       // MCMC settings.
		       unsigned int burnin, unsigned int interval, unsigned int dyninterval,
		       // Verbosity.
		       int fVerbose){


  unsigned int n_par=F_m->n_stats,
    n_stats=F_m->n_stats;

  double *F_stats_acc=malloc(sizeof(double)*n_stats),
    *F_stats2_acc=malloc(sizeof(double)*n_stats),
    *F_thetaP=malloc(sizeof(double)*n_par),
    *F_DobjDtheta=malloc(sizeof(double)*n_par),
    *delta=malloc(sizeof(double)*n_par),
    objPU,objPD,
    *D_stats=malloc(sizeof(double)*D_m->n_stats),
    *F_theta_acc=calloc(n_par,sizeof(double)),
    *F_theta2_acc=calloc(n_par,sizeof(double)),
    *F_theta_sd=calloc(n_par,sizeof(double)),
    acc_total=0,
    sdsum;

  int use_var=-20;

  // Burn-in
  MCMCSampleDynObjective(nwp, order, F_m, F_MH, F_theta, D_m, D_MH, D_theta, dev, D_stats, nmax, difftime, difftail, diffhead, dyninterval,burnin,interval,F_stats_acc,F_stats2_acc, &use_var, fVerbose);

  for(unsigned int i=0; i<iterations; i++){
    R_CheckUserInterrupt();
    double gain=a/pow(A+i+1,alpha);
    double diff=c/pow(i+1,gamma);
    
    if(fVerbose){
      Rprintf("\nIteration %d/%d: gain=%f diff=%f\n",i+1,iterations,gain,diff);
      Rprintf("F_theta=[ ");
      for(unsigned int k=0; k<n_par; k++) Rprintf("%f ",F_theta[k]);
      Rprintf("]\n");
    }

    sdsum=0;
    for(unsigned int k=0; k<n_par; k++){
      
      F_theta_acc[k]*=1-1.0/iterations*10;
      F_theta_acc[k]+=F_theta[k];
      F_theta2_acc[k]*=1-1.0/iterations*10;
      F_theta2_acc[k]+=F_theta[k]*F_theta[k];
      acc_total*=1-1.0/iterations*10;
      acc_total++;
      /*      if(i)
	F_theta_sd[k]=sqrt((F_theta2_acc[k]-F_theta_acc[k]*F_theta_acc[k]/(acc_total))/(acc_total));
	else*/
	F_theta_sd[k]=1;
      
      sdsum+=F_theta_sd[k];
    }

    // Generate delta
    if(fVerbose) Rprintf("Perturbation: [ ");
    for(unsigned int k=0; k<n_par; k++){
      delta[k]=(rbinom(1,0.5)*2-1)*diff*F_theta_sd[k]/sdsum*n_par;
      if(fVerbose) Rprintf("%f ", delta[k]);     
    }
    if(fVerbose) Rprintf("]\n");
    
    // Compute theta perturbed "up"
    for(unsigned int k=0; k<n_par; k++){
      F_thetaP[k]=F_theta[k]+delta[k];
    }
    // Evaluate the objective function with theta perturbed "up"
    objPU=MCMCSampleDynObjective(nwp, order, F_m, F_MH, F_thetaP, D_m, D_MH, D_theta, dev, D_stats, nmax, difftime, difftail, diffhead, dyninterval,interval,interval,F_stats_acc,F_stats2_acc, &use_var, fVerbose);

    if(use_var==0){
      if(fVerbose) Rprintf("Switching to variance-normalized objective function!\n");
      continue;
    }

    // Compute theta perturbed "down"
    for(unsigned int k=0; k<n_par; k++){
      F_thetaP[k]=F_theta[k]-delta[k];
    }
    // Evaluate the objective function with theta perturbed "down"
    objPD=MCMCSampleDynObjective(nwp, order, F_m, F_MH, F_thetaP, D_m, D_MH, D_theta, dev, D_stats, nmax, difftime, difftail, diffhead, dyninterval,interval,interval,F_stats_acc,F_stats2_acc, &use_var, fVerbose);

    if(use_var==0){
      if(fVerbose) Rprintf("Switching to variance-normalized objective function!\n");
    }

    // Estimate the gradient, and make the step (at most diff)
    if(fVerbose) Rprintf("Estimated gradient: [ ");
    for(unsigned int k=0; k<n_par; k++){
      F_DobjDtheta[k]=(objPU-objPD)/delta[k]/2;
      if(fVerbose) Rprintf("%f ", F_DobjDtheta[k]);
      F_theta[k]=F_theta[k]-sign(F_DobjDtheta[k])*fmin(fabs(gain*F_DobjDtheta[k]),diff);
      if(isnan(F_theta[k])) error("Dammit!");
    }
    if(fVerbose){
      Rprintf("]\n");
      Rprintf("F_theta=[ ");
      for(unsigned int k=0; k<n_par; k++) Rprintf("%f ",F_theta[k]);
      Rprintf("]\n");

      Rprintf("mean F_theta=[ ");
      for(unsigned int k=0; k<n_par; k++) Rprintf("%f ",F_theta_acc[k]/acc_total);
      Rprintf("]\n");
    }
    obj_history[i]=(objPU+objPD)/2;
  }

  // Return smoothed theta.
  for(unsigned int k=0; k<n_par; k++)
    F_theta[k]=F_theta_acc[k]/acc_total;

  free(F_stats_acc);
  free(F_stats2_acc);
  free(F_thetaP);
  free(F_DobjDtheta);
  free(F_theta_sd);
  free(D_stats);
  free(delta);
}


