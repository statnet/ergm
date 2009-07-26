#include "MCMC.h"
#include "MCMCDyn.h"
#include "MHproposals.h"

/*****************
 Note on undirected networks:  For j<k, edge {j,k} should be stored
 as (j,k) rather than (k,j).  In other words, only directed networks
 should have (k,j) with k>j.
*****************/

/*R_INLINE*/ void MCMCDyn_init_common(int *heads, int *tails, int n_edges,
          int maxedges,
				  int n_nodes, int dflag, int bipartite, Network *nw,

				  int order_code, DynamOrder *order,
				  int F_nterms, char *F_funnames, char *F_sonames, double *F_inputs, Model **F_m,
				  int D_nterms, char *D_funnames, char *D_sonames, double *D_inputs, Model **D_m,
				  
				  int *attribs, int *maxout, int *maxin, int *minout,
				  int *minin, int condAllDegExact, int attriblength, DegreeBound **bd,

				  char *F_MHproposaltype, char *F_MHproposalpackage, MHproposal *F_MH,
				  char *D_MHproposaltype, char *D_MHproposalpackage, MHproposal *D_MH,
				  int fVerbose){
  GetRNGstate();  /* R function enabling uniform RNG */
  
  switch(order_code){
  case 1: *order=DissThenForm; break;
  case 2: *order=DissAndForm; break;
  case 3: *order=FormThenDiss; break;
  case 4: *order=FormOnly; break;
  case 5: *order=DissOnly; break;
  default: error("Unsupported dynamic model code %d.", order_code);
  }
  if(fVerbose){
    switch(order_code){
    case 1: Rprintf("Using dissolve then form dynamic model code.\n"); break;
    case 2: Rprintf("Using simultaneous dissolve and form dynamic model code.\n"); break;
    case 3: Rprintf("Using form then dissolve dynamic model code.\n"); break;
    case 4: Rprintf("Using only formation dynamic model code.\n"); break;
    case 5: Rprintf("Using only dissolution dynamic model code.\n"); break;
    default: error("Unsupported dynamic model code %d.", order_code);
    }
  }

  *F_m=ModelInitialize(F_funnames, F_sonames, F_inputs, F_nterms);
  *D_m=ModelInitialize(D_funnames, D_sonames, D_inputs, D_nterms);

  nw[0]=NetworkInitialize(heads, tails, n_edges, 
                          n_nodes, dflag, bipartite, 1);
  nw[1]=NetworkInitialize(NULL, NULL, 0,
                          n_nodes, dflag, bipartite, 0);

  
  *bd=DegreeBoundInitialize(attribs, maxout, maxin, minout, minin,
			    condAllDegExact, attriblength, nw);

  MH_init(F_MH, F_MHproposaltype, F_MHproposalpackage, fVerbose, nw, *bd);
  MH_init(D_MH, D_MHproposaltype, D_MHproposalpackage, fVerbose, nw, *bd);

}

/*R_INLINE*/ void MCMCDyn_finish_common(Network *nw,

				    Model *F_m,
				    Model *D_m,
				  
				    DegreeBound *bd,

				    MHproposal *F_MH,
				    MHproposal *D_MH){
  MH_free(F_MH);
  MH_free(D_MH);
  if(bd)DegreeBoundDestroy(bd);
  ModelDestroy(F_m);
  ModelDestroy(D_m);
  NetworkDestroy(nw);
  NetworkDestroy(&nw[1]);
  PutRNGstate();  /* Disable RNG before returning */

}

/*****************
 void MCMC_wrapper

 Wrapper for a call from R.
*****************/
void MCMCDyn_wrapper(// Starting network.
		     int *heads, int *tails, int *n_edges, int *maxpossibleedges,
		     int *n_nodes, int *dflag, int *bipartite,
		     // Ordering of formation and dissolution.
		     int *order_code,
		     // Formation terms and proposals.
		     int *F_nterms, char **F_funnames, char **F_sonames, 
		     char **F_MHproposaltype, char **F_MHproposalpackage,
		     double *F_inputs, double *theta, 
		     // Dissolution terms and proposals.
		     int *D_nterms, char **D_funnames, char **D_sonames,
		     char **D_MHproposaltype, char **D_MHproposalpackage,
		     double *D_inputs, double *gamma,
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
		     int *fVerbose){
  int i;
  Edge  nmax;
  Network nw[2];
  DegreeBound *bd;
  Model *F_m, *D_m;
  MHproposal F_MH, D_MH;
  DynamOrder order;

  nmax = (Edge)*maxedges; /* coerce double *maxedges to type Edge */
  
  
  Vertex *difftime, *diffhead, *difftail;
  difftime = (Vertex *) diffnetworktime;
  diffhead = (Vertex *) diffnetworkhead;
  difftail = (Vertex *) diffnetworktail;

  for (i = 0; i < nmax; i++){
    newnetworkhead[i] = 0;
    newnetworktail[i] = 0;
    if(*burnin==0 && *interval==1){
      difftime[i] = 0;
      diffhead[i] = 0;
      difftail[i] = 0;
    }
  }

  MCMCDyn_init_common(heads, tails, *n_edges, *maxpossibleedges,
		      *n_nodes, *dflag, *bipartite, nw,
		      *order_code, &order,
		      *F_nterms, *F_funnames, *F_sonames, F_inputs, &F_m,
		      *D_nterms, *D_funnames, *D_sonames, D_inputs, &D_m,
		      attribs, maxout, maxin, minout,
		      minin, *condAllDegExact, *attriblength, &bd,
		      *F_MHproposaltype, *F_MHproposalpackage, &F_MH,
		      *D_MHproposaltype, *D_MHproposalpackage, &D_MH,
		      *fVerbose);

  MCMCSampleDyn(nw, order,
		F_m, &F_MH, theta,
		D_m, &D_MH, gamma,
		bd,
		F_sample, D_sample, nmax, difftime, diffhead, difftail,
		*nsteps, *dyninterval, *burnin, *interval,
		*fVerbose);
   
  /* record new generated network to pass back to R */

  newnetworktail[0]=newnetworkhead[0]=EdgeTree2EdgeList(newnetworkhead+1,newnetworktail+1,nw,nmax);

  MCMCDyn_finish_common(nw, F_m, D_m, bd, &F_MH, &D_MH);

}

/*********************
 void MCMCSampleDyn

 Using the parameters contained in the array theta, obtain the
 network statistics for a sample of size nsteps.  burnin is the
 initial number of Markov chain steps before sampling anything
 and interval is the number of MC steps between successive 
 networks in the sample.  Put all the sampled statistics into
 the F_statistics array. 
*********************/
void MCMCSampleDyn(// Observed and discordant network.
		   Network *nwp,
		   // Ordering of formation and dissolution.
		   DynamOrder order,
		   // Formation terms and proposals.
		   Model *F_m, MHproposal *F_MH,
		   double *theta,
		   // Dissolution terms and proposals.
		   Model *D_m, MHproposal *D_MH,
		   double *gamma,
		   // Degree bounds.
		   DegreeBound *bd,
		   // Space for output.
		   double *F_stats, double *D_stats, // Do we still need these?
		   Edge nmax,
		   Vertex *difftime, Vertex *diffhead, Vertex *difftail,		    
		   // MCMC settings.
		   unsigned int nsteps, unsigned int dyninterval,
		   unsigned int burnin, unsigned int interval, 
		   // Verbosity.
		   int fVerbose){

  int i, j, log_toggles = (burnin==0 && interval==1);
  Edge nextdiffedge=1;

  if (fVerbose)
    Rprintf("Total m->n_stats is %i; total nsteps is %d\n",
	    F_m->n_stats,nsteps);
  
  
  /* Burn in step. */

  for(i=0;i<burnin;i++)
    MCMCDyn1Step(nwp, order, 
		 F_m, F_MH, theta, D_m, D_MH, gamma, bd,
		 log_toggles, F_stats, D_stats,
		 nmax, &nextdiffedge, difftime, diffhead, difftail,
		 dyninterval, fVerbose);
  
  //Rprintf("MCMCSampleDyn post burnin numdissolve %d\n", *numdissolve);
  
  if (fVerbose){
    Rprintf("Returned from Metropolis-Hastings burnin\n");
  }
  
  /* Now sample networks */
  for (i=0; i < nsteps; i++){
    /* Set current vector of stats equal to previous vector */
    for (j=0; j<F_m->n_stats; j++){
      F_stats[j+F_m->n_stats] = F_stats[j];
    }
    F_stats += F_m->n_stats;

    for (j=0; j<D_m->n_stats; j++){
      D_stats[j+D_m->n_stats] = D_stats[j];
    }
    D_stats += D_m->n_stats;

    /* This then adds the change statistics to these values */
    for(j=0;j<interval;j++){
      MCMCDyn1Step(nwp, order, 
		   F_m, F_MH, theta, D_m, D_MH, gamma, bd,
		   log_toggles, F_stats, D_stats,
		   nmax, &nextdiffedge, difftime, diffhead, difftail,
		   dyninterval, fVerbose);
      if(log_toggles && nextdiffedge>=nmax) {
	if(fVerbose) Rprintf("Nmax of %d exceeded.\n",nmax);
	// Mark this run as "bad".
	difftime[0]=diffhead[0]=difftail[0]=nmax+1;
	return;
      }
    }
    
    //Rprintf("MCMCSampleDyn loop numdissolve %d\n", *numdissolve);
    if (fVerbose){
      if( ((3*i) % nsteps)<3 && nsteps > 500){
	Rprintf("Advanced %d time steps.\n", i);
      }
    }
  }

  if(log_toggles) difftime[0]=diffhead[0]=difftail[0]=nextdiffedge-1;
}

/* Helper function to run the formation or the dissolution side of the process, 
   depending on proposal, statistics, and parameters passed. 

   NOTE: Here, "stats" and "m" are for the process that "drives" the sampling 
   in this phase, while "O_stats" and "O_m" are for the other process (whose
   statistcs are still affected and need to be updated).
*/
/*R_INLINE*/ void MCMCDyn1Step_sample(MHproposal *MH,
				  double *par, double *stats,
				  int dyninterval, 
				  Network *nwp,
				  Model *m, DegreeBound *bd){
  Vertex step;
  double cutoff, ip;
  unsigned int i;
  
  MH->ntoggles = 0;
  (*(MH->func))(MH, bd, nwp); /* Call MH proposal function to initialize */
  
  for(step = 0; step < dyninterval; step++) {
    
    MH->ratio = 1.0;
    (*(MH->func))(MH, bd, nwp); /* Call MH function to propose toggles */
    //      Rprintf("Back from proposal; step=%d\n",step);

    // Proposal failed.
    if(MH->togglehead[0]==MH_FAILED){
      if(MH->toggletail[0]==MH_UNRECOVERABLE) 
	error("Something very bad happened during proposal. Memory has not been deallocated, so restart R soon.");
      if(MH->toggletail[0]==MH_IMPOSSIBLE) break;
    }

    ChangeStats(MH->ntoggles, MH->togglehead, MH->toggletail, nwp, m);
    
    //  Rprintf("change stats:"); 
    /* Calculate inner product */
    for (i=0, ip=0.0; i<m->n_stats; i++){
      ip += par[i] * m->workspace[i];
      //  Rprintf("%f ", m->workspace[i]); 
    }
    //  Rprintf("\n ip %f dedges %f\n", ip, m->workspace[0]); 
    /* The logic is to set exp(cutoff) = exp(ip) * qratio ,
       then let the MH probability equal min{exp(cutoff), 1.0}.
       But we'll do it in log space instead.  */
    cutoff = ip + log(MH->ratio);
    
    /* if we accept the proposed network */
    if (cutoff >= 0.0 || log(unif_rand()) < cutoff) { 
      /* Hold off updating timesteamps until the changes are committed,
	 which doesn't happen until later. */
      for (i=0; i < MH->ntoggles; i++){
	ToggleEdge(MH->togglehead[i], MH->toggletail[i], &nwp[0]);
	ToggleEdge(MH->togglehead[i], MH->toggletail[i], &nwp[1]);  /* Toggle the discord for this edge */
      }
      /* Do NOT record network statistics for posterity yet. */
    }
  }
}

/*
  MCMCDyn1Step_commit
  Applies a list of toggles to a network, updating change statistics and time stamps.
*/
/*R_INLINE*/ void MCMCDyn1Step_commit(unsigned int ntoggles,
				  Vertex *diffhead, Vertex *difftail,
				  Network *nwp,
				  Model *F_m, double *F_stats,
				  Model *D_m, double *D_stats){
  ChangeStats(ntoggles,diffhead,difftail,nwp,F_m);
  ChangeStats(ntoggles,diffhead,difftail,nwp,D_m);
  for (unsigned int i = 0; i < F_m->n_stats; i++)
    F_stats[i] += F_m->workspace[i];
  for (unsigned int i = 0; i < D_m->n_stats; i++)
    D_stats[i] += D_m->workspace[i];

  for(Edge i=0;i<ntoggles;i++){
    ToggleEdgeWithTimestamp(diffhead[i],difftail[i],nwp);
  }
}

/* 
   MCMCDyn1Step_record_reset
   Record a set of changes to the state of the sampler at the end of a phase:
   * record new generated network differences to pass back to R
   * undo the toggles to nwp[0]
   * empty the discordant network (nwp[1])
   * return the number of edges toggled
*/
/*R_INLINE*/ unsigned int MCMCDyn1Step_record_reset(Edge nmax,
						Vertex *difftime, Vertex *diffhead, Vertex *difftail,
						Network *nwp, 
						Edge *nextdiffedge){
  Vertex head, tail;
  const unsigned int t=nwp->duration_info.MCMCtimer;
  Edge ntoggles = nwp[1].nedges;
  
  for(unsigned int i=0; i<ntoggles; i++){
    FindithEdge(&head, &tail, 1, &nwp[1]); // Grab the next edge that changed;
    ToggleEdge(head, tail, &nwp[1]); // delete it from nwp[1];
    ToggleEdge(head, tail, nwp); // undo the toggle in nwp[0];
    
    if(*nextdiffedge<nmax){
      // and record the toggle.
      difftime[*nextdiffedge] = t;
      diffhead[*nextdiffedge] = head;
      difftail[*nextdiffedge] = tail;
      (*nextdiffedge)++;
    }
  }
  return(ntoggles);
}

/*********************
 void MCMCDyn1Step

 Simulate evolution of a dynamic network for nsteps steps.
*********************/
void MCMCDyn1Step(// Observed and discordant network.
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
		  unsigned log_toggles,
		  double *F_stats, double *D_stats,
		  unsigned int nmax, Edge *nextdiffedge,
		  Vertex *difftime, Vertex *diffhead, Vertex *difftail,
		  // MCMC settings.
		  unsigned int dyninterval,
		  // Verbosity.
		  int fVerbose){
  
  Edge ntoggles;
  Edge startnextdiffedge=*nextdiffedge;

  /* Increment the MCMC timer. */
  nwp->duration_info.MCMCtimer++;

  switch(order){
  case DissThenForm:
    /* Run the dissolution process and commit it. */
    MCMCDyn1Step_sample(D_MH, gamma, D_stats, dyninterval, nwp, D_m, bd);
    ntoggles = MCMCDyn1Step_record_reset(nmax, difftime, diffhead, difftail, nwp, nextdiffedge);
    MCMCDyn1Step_commit(ntoggles, diffhead+*nextdiffedge-ntoggles, difftail+*nextdiffedge-ntoggles, nwp, F_m, F_stats, D_m, D_stats);
    
    /* Run the formation process and commit it. */
    MCMCDyn1Step_sample(F_MH, theta, F_stats, dyninterval, nwp, F_m, bd);
    ntoggles = MCMCDyn1Step_record_reset(nmax, difftime, diffhead, difftail, nwp, nextdiffedge);
    MCMCDyn1Step_commit(ntoggles, diffhead+*nextdiffedge-ntoggles, difftail+*nextdiffedge-ntoggles, nwp, F_m, F_stats, D_m, D_stats);
    
    break;
  case DissAndForm:
    /* Run the dissolution process. */
    MCMCDyn1Step_sample(D_MH, gamma, D_stats, dyninterval, nwp, D_m, bd);
    ntoggles = MCMCDyn1Step_record_reset(nmax, difftime, diffhead, difftail, nwp, nextdiffedge);
    
    /* Run the formation process. */
    MCMCDyn1Step_sample(F_MH, theta, F_stats, dyninterval, nwp, F_m, bd);
    ntoggles += MCMCDyn1Step_record_reset(nmax, difftime, diffhead, difftail, nwp, nextdiffedge);
    
    /* Commit both. */
    MCMCDyn1Step_commit(ntoggles, diffhead+*nextdiffedge-ntoggles, difftail+*nextdiffedge-ntoggles, nwp, F_m, F_stats, D_m, D_stats);
    break;
  case FormThenDiss:
    /* Run the formation process and commit it. */
    MCMCDyn1Step_sample(F_MH, theta, F_stats, dyninterval, nwp, F_m, bd);
    ntoggles = MCMCDyn1Step_record_reset(nmax, difftime, diffhead, difftail, nwp, nextdiffedge);
    MCMCDyn1Step_commit(ntoggles, diffhead+*nextdiffedge-ntoggles, difftail+*nextdiffedge-ntoggles, nwp, F_m, F_stats, D_m, D_stats);

    /* Run the dissolution process and commit it. */
    MCMCDyn1Step_sample(D_MH, gamma, D_stats, dyninterval, nwp, D_m, bd);
    ntoggles = MCMCDyn1Step_record_reset(nmax, difftime, diffhead, difftail, nwp, nextdiffedge);
    MCMCDyn1Step_commit(ntoggles, diffhead+*nextdiffedge-ntoggles, difftail+*nextdiffedge-ntoggles, nwp, F_m, F_stats, D_m, D_stats);
    break;
  case FormOnly:
    /* Run the formation process and commit it. */
    MCMCDyn1Step_sample(F_MH, theta, F_stats, dyninterval, nwp, F_m, bd);
    ntoggles = MCMCDyn1Step_record_reset(nmax, difftime, diffhead, difftail, nwp, nextdiffedge);
    MCMCDyn1Step_commit(ntoggles, diffhead+*nextdiffedge-ntoggles, difftail+*nextdiffedge-ntoggles, nwp, F_m, F_stats, D_m, D_stats);
    break;
  case DissOnly:
    /* Run the dissolution process and commit it. */
    MCMCDyn1Step_sample(D_MH, gamma, D_stats, dyninterval, nwp, D_m, bd);
    ntoggles = MCMCDyn1Step_record_reset(nmax, difftime, diffhead, difftail, nwp, nextdiffedge);
    MCMCDyn1Step_commit(ntoggles, diffhead+*nextdiffedge-ntoggles, difftail+*nextdiffedge-ntoggles, nwp, F_m, F_stats, D_m, D_stats);
    break;
  default: 
    error("Unsupported dynamic model code %d. Memory has not been deallocated so restart R sometime soon.",order);
  }

  // If we don't keep a log of toggles, reset the position to save space.
  if(!log_toggles)
    *nextdiffedge=startnextdiffedge;
}

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
		      *D_nterms, *D_funnames, *D_sonames, D_inputs, &D_m,
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
  unsigned int phase1n=phase1n_base+3*F_m->n_stats, *changed;
  
  double *ubar, *aDdiaginv, *D_stats, *prevdev, n, ubar0, ubar0_next;

  unsigned int nomix;
  unsigned int nomixed=0;
  
  do{
    nomixed++;
    if(nomixed>10) error("Robbins-Monro failing to converge.");
    ubar = (double *)malloc( F_m->n_stats * sizeof(double));
    prevdev = (double *)malloc( F_m->n_stats * sizeof(double));
    changed = (unsigned int *)malloc( F_m->n_stats * sizeof(double));
    aDdiaginv = (double *)malloc( F_m->n_stats * sizeof(double));
    D_stats = (double *)calloc( D_m->n_stats, sizeof(double));
    
    Edge nextdiffedge=1;

    ubar0_next=0;

    for (j=0; j < F_m->n_stats; j++){
      ubar[j] = 0.0;
      prevdev[j]=dev[j];
      changed[j]=0;
      Rprintf("j %d %f\n",j,theta[j]);
      n=0;
      theta[j]-=0.5;
    }
    
    /*********************
    Burn in step. 
    *********************/
    
    if(fVerbose) Rprintf("Starting burnin of %d steps\n", burnin);
    for(i=0;i<burnin;i++){
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

      // Find the average value for the first statistic.
      ubar0_next *= 1-1.0/burnin;
      n*=1-1.0/burnin;
      ubar0_next += dev[0];
      n++;
      
      for(j=0; j<F_m->n_stats; j++){
	if(i>burnin/2 && dev[j]!=prevdev[j]) changed[j]++;
	prevdev[j]=dev[j];
      }
    }

    ubar0_next /= n;
    nomix=0;
    for(j=0; j<F_m->n_stats; j++){
            if(changed[j]<0.05*burnin/2 && dev[j]!=0){ /* It's OK if statistic is spot on. */
	if(fVerbose)Rprintf("Bad mixing: %d\n", j);
	nomix=1;
      }
    }
    
    /********************
      Phase 1: estimate dgy/dtheta
    ********************/
    Rprintf("Phase 1: %d steps (interval = %d)\n", phase1n,interval);
    
    for(j=0; j<F_m->n_stats; j++){
      ubar0=ubar0_next;
      ubar0_next=0;

      /*If it's a badly mixing statistic, don't bother with its derivatie this round.*/
      if(changed[j]<0.05*burnin/2 && dev[j]!=0){
	ubar[j]=0;
	continue;
      }

      theta[j]++;
      n=0;
      for(i=0; i < phase1n*interval; i++){
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
	ubar[j]*=1-1.0/(phase1n*interval);
	n*=1-1.0/(phase1n*interval);
	ubar[j] += dev[j]-ubar0;
	n++;
	
	if(j+1<F_m->n_stats){
	  ubar0_next*=1-1.0/(phase1n*interval);
	  ubar0_next+=dev[j+1];
	}
      }
      
      ubar[j]/=n;
      ubar0_next/=n;
    }
    
    if (fVerbose){
      Rprintf("Returned from Phase 1\n");
      Rprintf("j, approx dgy/dtheta, gain*dtheta/dgy:\n");
    }
    
    for (j=0; j<F_m->n_stats; j++){
      aDdiaginv[j]=ubar[j];
      if( aDdiaginv[j] > 0.0){
	if(fVerbose) Rprintf("%d, %f, %f\n", j, aDdiaginv[j], gain/aDdiaginv[j]);
	aDdiaginv[j] = gain/aDdiaginv[j];
      }else{
	if(fVerbose) Rprintf("%d, %f, %f\n", j, aDdiaginv[j], 0.00000);
	aDdiaginv[j]=0.00000;
	nomix=1;
      }
    }
    
    /********************
      Phase 2
    ********************/
    unsigned int phase2n=F_m->n_stats+7+phase2n_base;
    double *meandev=(double*)calloc(F_m->n_stats,sizeof(double));
    double *meandevlong=(double*)calloc(F_m->n_stats,sizeof(double));
    double n2;
    
    for(unsigned int subphase=0; subphase<(nomix?(int)ceil(phase2sub/2):phase2sub); subphase++){
      
      for(j=0; j<F_m->n_stats; j++) meandevlong[j]=dev[j];
      n2=1;
      for (i=0; i < phase2n; i++){
	for(j=0; j<F_m->n_stats; j++) meandev[j]=dev[j];
	n=1;
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
	  for(unsigned int k=0;k<F_m->n_stats; k++){
	    meandev[k]*=1-1.0/interval;
	    n*=1-1.0/interval;
	    meandev[k]+=dev[k];
	    n++;
	  }
	}
	
	/* Update theta0 */
	for (j=0; j<F_m->n_stats; j++){
	  theta[j] -= aDdiaginv[j] * (meandev[j]/n);
	  
	  meandevlong[j]*=1-1.0/phase2n;
	  n2*=1-1.0/phase2n;
	  meandevlong[j]+=meandev[j]/n;
	  n2++;
	}
      }
      
      for (j=0; j<F_m->n_stats; j++){
	aDdiaginv[j] /= 2.0;
	if (fVerbose)Rprintf("subphase j %d theta %f ns %f\n",
			     j, theta[j], meandevlong[j]/n2);
      }
      Rprintf("\n");
      
      phase2n=trunc(2.52*(phase2n-phase2n_base)+phase2n_base);
    }
    
    
    free(ubar);
    free(prevdev);
    free(changed);
    free(aDdiaginv);
    free(D_stats );
  }while(nomix);
	 

}
