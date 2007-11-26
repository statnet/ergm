#include "MCMC.h"
#include "MCMCDyn.h"
#include "MHproposals.h"

/*****************
 Note on undirected networks:  For j<k, edge {j,k} should be stored
 as (j,k) rather than (k,j).  In other words, only directed networks
 should have (k,j) with k>j.
*****************/

/*****************
 void MCMC_wrapper

 Wrapper for a call from R.
*****************/
void MCMCDyn_wrapper(// Starting network.
		     int *heads, int *tails, int *n_edges,
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

  if((diffnetworktime || diffnetworkhead || diffnetworktail) &&
     (*burnin!=0 || *interval!=1)){
    if(fVerbose) Rprintf("WARNING: Burnin is not 0 or interval is not 1. Diff lists are not meaningful.\n");
  }

  switch(*order_code){
  case 1: order=DissThenForm; break;
  case 2: order=DissAndForm; break;
  case 3: order=FormThenDiss; break;
  case 4: order=FormOnly; break;
  case 5: order=DissOnly; break;
  default: error("Unsupported dynamic model code %d.", *order_code);
  }
  if(*fVerbose){
    switch(*order_code){
    case 1: Rprintf("Using dissolve then form dynamic model code.\n"); break;
    case 2: Rprintf("Using simultaneous dissolve and form dynamic model code.\n"); break;
    case 3: Rprintf("Using form then dissolve dynamic model code.\n"); break;
    case 4: Rprintf("Using only formation dynamic model code.\n"); break;
    case 5: Rprintf("Using only dissolution dynamic model code.\n"); break;
    default: error("Unsupported dynamic model code %d.", *order_code);
    }
  }
  
  nmax = (Edge)*maxedges; /* coerce double *maxedges to type Edge */
  
  GetRNGstate();  /* R function enabling uniform RNG */
  
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

  F_m=ModelInitialize(*F_funnames, *F_sonames, F_inputs, *F_nterms);
  D_m=ModelInitialize(*D_funnames, *D_sonames, D_inputs, *D_nterms);

  nw[0]=NetworkInitialize(heads, tails, *n_edges, *n_nodes, *dflag, *bipartite, 1);
  nw[1]=NetworkInitialize(NULL, NULL, 0, *n_nodes, *dflag, *bipartite, 0);


  bd=DegreeBoundInitialize(attribs, maxout, maxin, minout, minin,
			   *condAllDegExact, *attriblength, nw);

  MH_init(&F_MH, *F_MHproposaltype, *F_MHproposalpackage, *fVerbose, nw, bd);
  MH_init(&D_MH, *D_MHproposaltype, *D_MHproposalpackage, *fVerbose, nw, bd);

  MCMCSampleDyn(nw, order,
		F_m, &F_MH, theta,
		D_m, &D_MH, gamma,
		bd,
		F_sample, D_sample, nmax, difftime, diffhead, difftail,
		*nsteps, *dyninterval, *burnin, *interval,
		*fVerbose);
   
// Rprintf("nsteps: %f\n", *nsteps);
// for (i=0; i < (*nsteps)*9; i++){
// 	if(i == 9*trunc(i/9)){Rprintf("\n");}
// Rprintf("%f ", sample[i]);
// }
// Rprintf("\n");

  /* record new generated network to pass back to R */

  newnetworktail[0]=newnetworkhead[0]=EdgeTree2EdgeList(newnetworkhead+1,newnetworktail+1,nw,nmax);

  MH_free(&F_MH);
  MH_free(&D_MH);
  DegreeBoundDestroy(bd);
  ModelDestroy(F_m);
  ModelDestroy(D_m);
  NetworkDestroy(nw);
  NetworkDestroy(&nw[1]);

  PutRNGstate();  /* Disable RNG before returning */
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

  int i, j;
  Edge nextdiffedge=1;

  if (fVerbose)
    Rprintf("Total m->n_stats is %i; total nsteps is %d\n",
	    F_m->n_stats,nsteps);
  
  
  /* Burn in step. */

  for(i=0;i<burnin;i++)
    MCMCDyn1Step(nwp, order, 
		 F_m, F_MH, theta, D_m, D_MH, gamma, bd,
		 0, F_stats, D_stats,
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
		   burnin==0 && interval==1, F_stats, D_stats,
		   nmax, &nextdiffedge, difftime, diffhead, difftail,
		   dyninterval, fVerbose);
      if(nextdiffedge>=nmax) {
	if(fVerbose) Rprintf("Nmax of %d exceeded.\n",nmax);
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

  if(burnin==0 && interval==1) difftime[0]=diffhead[0]=difftail[0]=nextdiffedge-1;
}

/* Helper function to run the formation or the dissolution side of the process, 
   depending on proposal, statistics, and parameters passed. 

   NOTE: Here, "stats" and "m" are for the process that "drives" the sampling 
   in this phase, while "O_stats" and "O_m" are for the other process (whose
   statistcs are still affected and need to be updated).
*/
R_INLINE void MCMCDyn1Step_sample(MHproposal *MH,
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
      if(MH->togglehead[1]==MH_UNRECOVERABLE) 
	error("Something very bad happened during proposal. Memory has not been deallocated, so restart R soon.");
      if(MH->togglehead[1]==MH_IMPOSSIBLE) break;
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
R_INLINE void MCMCDyn1Step_commit(unsigned int ntoggles,
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
R_INLINE unsigned int MCMCDyn1Step_record_reset(Edge nmax,
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
  
  //FIXME: log_toggles directive is ignored!

  Edge ntoggles;

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
}

void MCMCDynPhase12(// Observed network.
		    int *heads, int *tails, int *n_edges,
		    int *n_nodes, int *dflag, int *bipartite, 
		    // Ordering of formation and dissolution.
		    int *order_code, 
		    // Formation terms and proposals.
		    int *F_nterms, char **F_funnames, char **F_sonames, 
		    char **F_MHproposaltype, char **F_MHproposalpackage,
		    double *F_inputs, double *theta0, 
		    // Formation parameter fitting.
		    double *gain, double *meanstats, 
		    int *phase1, int *nsub,
		    // Dissolution terms and proposals.
		    int *D_nterms, char **D_funnames, char **D_sonames, 
		    char **D_MHproposaltype, char **D_MHproposalpackage,
		    double *D_inputs, double *gamma0,
		    // Dissolution parameter fitting --- to add later? -PK
		    // Degree bounds.
		    int *attribs, int *maxout, int *maxin, int *minout,
		    int *minin, int *condAllDegExact, int *attriblength, 
		    // MCMC settings.
		    double *nsteps, int *dyninterval,
		    double *burnin, double *interval, 
		    // Space for output.
		    double *F_sample, double *D_sample, 
		    int *newnetworkhead, int *newnetworktail, 
		    double *maxedges,
		    int *diffnetworktime, int *diffnetworkhead, int *diffnetworktail, 
		    // Verbosity.
		    int *fVerbose){
  int i, nextedge, directed_flag, hammingterm, formationterm;
  int nphase1, nsubphases;
  Vertex v, k, bip, hhead, htail;
  Edge nddyads, kedge, nmax;
  Network nw[2];
  DegreeBound *bd;
  Model *F_m, *D_m;
  MHproposal F_MH, D_MH;
  ModelTerm *thisterm;
  DynamOrder order;

  switch(*order_code){
  case 1: order=DissThenForm; break;
  case 2: order=DissAndForm; break;
  case 3: order=FormThenDiss; break;
  default: error("Unsupported dynamic model code %d.", order_code);
  }
  
  nphase1 = (int)*phase1; /* coerce double *n_nodes to type Vertex */
  nsubphases = (int)*nsub; /* coerce double *n_nodes to type Vertex */
  nmax = (Edge)*maxedges; /* coerce double *maxedges to type Edge */
  bip = (Vertex)*bipartite; /* coerce double *bipartite to type Vertex */
  
  GetRNGstate();  /* R function enabling uniform RNG */
  
  directed_flag = *dflag;

  Vertex *difftime, *diffhead, *difftail;
  difftime = (Vertex *) diffnetworktime;
  diffhead = (Vertex *) diffnetworkhead;
  difftail = (Vertex *) diffnetworktail;

  for (i = 0; i < nmax; i++){
    newnetworkhead[i] = 0;
    newnetworktail[i] = 0;
    difftime[i] = 0;
    diffhead[i] = 0;
    difftail[i] = 0;
  }

  F_m=ModelInitialize(*F_funnames, *F_sonames, F_inputs, *F_nterms);
  D_m=ModelInitialize(*D_funnames, *D_sonames, D_inputs, *D_nterms);

  nw[0]=NetworkInitialize(heads, tails, *n_edges, *n_nodes, *dflag, *bipartite, 1);
  nw[1]=NetworkInitialize(NULL, NULL, 0, *n_nodes, *dflag, *bipartite, 0);

  bd=DegreeBoundInitialize(attribs, maxout, maxin, minout, minin,
			   *condAllDegExact, *attriblength, nw);

  MH_init(&F_MH, *F_MHproposaltype, *F_MHproposalpackage, *fVerbose, nw, bd);
  MH_init(&D_MH, *D_MHproposaltype, *D_MHproposalpackage, *fVerbose, nw, bd);
  
//Rprintf("nsubphases %d\n", nsubphases);

  MCMCSampleDynPhase12(nw,order,
		       F_m, &F_MH, theta0, *gain, meanstats, nphase1, nsubphases, 
		       D_m, &D_MH, gamma0,
		       bd,
		       F_sample, D_sample,
		       nmax,
		       difftime, diffhead, difftail,
		       *nsteps, *burnin, *interval, *dyninterval,
		       *fVerbose);
   
// Rprintf("nsteps: %f\n", *nsteps);
// for (i=0; i < (*nsteps)*9; i++){
// 	if(i == 9*trunc(i/9)){Rprintf("\n");}
// Rprintf("%f ", sample[i]);
// }
// Rprintf("\n");

  /* record new generated network to pass back to R */

  newnetworktail[0]=newnetworkhead[0]=EdgeTree2EdgeList(newnetworkhead+1,newnetworktail+1,nw,nmax);

  MH_free(&F_MH);
  MH_free(&D_MH);
  DegreeBoundDestroy(bd);
  ModelDestroy(F_m);
  ModelDestroy(D_m);
  NetworkDestroy(nw);
  NetworkDestroy(&nw[1]);
  PutRNGstate();  /* Disable RNG before returning */
}


/*********************
 void MCMCSampleDynPhase12

 Using the parameters contained in the array theta, obtain the
 network statistics for a sample of size samplesize.  burnin is the
 initial number of Markov chain steps before sampling anything
 and interval is the number of MC steps between successive 
 networks in the sample.  Put all the sampled statistics into
 the F_stats array. 
*********************/
void MCMCSampleDynPhase12(// Observed and discordant network.
			  Network *nwp,
			  // Ordering of formation and dissolution.
			  DynamOrder order,
			  // Formation terms and proposals.
			  Model *F_m, MHproposal *F_MH,
			  double *theta, 
			  // Formation parameter fitting.
			  double gain, double *meanstats, int nphase1, int nsubphases,
			  // Dissolution terms and proposals.
			  Model *D_m, MHproposal *D_MH,
			  double *gamma, 
			  // Dissolution parameter fitting --- to add later? -PK
			  // Degree bounds.
			  DegreeBound *bd,
			  // Space for output.
			  double *F_stats,  double *D_stats, // Do we still need theese?
			  Edge nmax,
			  Vertex *difftime, Vertex *diffhead, Vertex *difftail,
			  // MCMC settings.
			  unsigned int samplesize, unsigned int burnin, 
			  unsigned int interval, unsigned int dyninterval,
			  // Verbosity.
			  int fVerbose){
  int i, j, components, diam;
  Edge nextdiffedge=1;
  ModelTerm *mtp;
  
  
//Rprintf("nsubphases %d\n", nsubphases);

  components = diam = 0;
  
  if (fVerbose)
    Rprintf("Total m->n_stats is %i; total samplesize is %d\n",
             F_m->n_stats,samplesize);

  /*********************
  F_stats are modified in groups of m->n_stats, and they
  reflect the CHANGE in the values of the statistics from the
  original (observed) network.  Thus, when we begin, the initial 
  values of the first group of m->n_stats F_stats should 
  all be zero
  *********************/
  double *ubar, *u2bar, *aDdiaginv;
  ubar = (double *)malloc( F_m->n_stats * sizeof(double));
  u2bar = (double *)malloc( F_m->n_stats * sizeof(double));
  aDdiaginv = (double *)malloc( F_m->n_stats * sizeof(double));
  Rprintf("\n");
  for (j=0; j < F_m->n_stats; j++){
    F_stats[j] = -meanstats[j];
    ubar[j] = 0.0;
    u2bar[j] = 0.0;
    Rprintf("j %d %f\n",j,theta[j]);
  }
  mtp = F_m->termarray;

  /*********************
   Burn in step.  While we're at it, use burnin statistics to 
   prepare covariance matrix for Mahalanobis distance calculations 
   in subsequent calls to M-H
   *********************/
//Rprintf("MCMCSampleDyn pre burnin numdissolve %d\n", *numdissolve);
  
    Rprintf("Starting burnin of %d steps\n", burnin);
    for(i=0;i<burnin;i++)
      MCMCDyn1Step(nwp, order,
		   F_m, F_MH, theta,
		   D_m, D_MH, gamma,
		   bd,
		   1,
		   F_stats, D_stats,
		   nmax, &nextdiffedge,
		   difftime, diffhead, difftail,
		   dyninterval,
		   fVerbose);

    Rprintf("Phase 1: %d steps (interval = %d)\n", nphase1,interval);
    /* Now sample networks */
    for (i=0; i <= nphase1; i++){
      for(j=0;j<interval;j++)
	MCMCDyn1Step(nwp, order,
		     F_m, F_MH, theta,
		     D_m, D_MH, gamma,
		     bd,
		     1,
		     F_stats, D_stats,
		     nmax, &nextdiffedge,
		     difftime, diffhead, difftail,
		     dyninterval,
		     fVerbose);
      if(i > 0){
	for (j=0; j<F_m->n_stats; j++){
	  ubar[j]  += F_stats[j];
	  u2bar[j] += F_stats[j]*F_stats[j];
	}
      }
      //  Rprintf("done %d step gain %f \n", i, gain);
    }
    if (fVerbose){
      Rprintf("Returned from Phase 1\n");
    }
    Rprintf("\n gain times inverse variances:\n");
    for (j=0; j<F_m->n_stats; j++){
      aDdiaginv[j] = u2bar[j]-ubar[j]*ubar[j]/(1.0*nphase1);
      if( aDdiaginv[j] > 0.0){
        aDdiaginv[j] = nphase1*gain/aDdiaginv[j];
      }else{
	aDdiaginv[j]=0.00001;
      }
      Rprintf(" %f", aDdiaginv[j]);
    }
    Rprintf("\n");
    
    /* Now sample networks */
    for (i=1; i < samplesize; i++){
      for(j=0;j<interval;j++)
	MCMCDyn1Step(nwp, order,
		     F_m, F_MH, theta,
		     D_m, D_MH, gamma,
		     bd,
		     1,
		     F_stats, D_stats,
		     nmax, &nextdiffedge,
		     difftime, diffhead, difftail,
		     dyninterval,
		     fVerbose);
      /* Update theta0 */
      for (j=0; j<F_m->n_stats; j++){
        theta[j] -= aDdiaginv[j] * F_stats[j];
      }
//Rprintf("\n");
//    if (fVerbose){ Rprintf("nsubphases %d i %d\n", nsubphases, i); }
      if (i==nsubphases){
	nsubphases = trunc(nsubphases*2.52) + 1;
        if (fVerbose){Rprintf("Updating nsub to be %d\n",nsubphases);}
        for (j=0; j<F_m->n_stats; j++){
          aDdiaginv[j] /= 2.0;
          if (fVerbose){Rprintf("j %d theta %f ns %f\n",
				j, theta[j], F_stats[j]);}
//        if (fVerbose){ Rprintf(" %f statsmean %f",  theta[j],(F_stats[j]-meanstats[j])); }
        }
        Rprintf("\n");
      }
      /* Set current vector of stats equal to previous vector */
      for (j=0; j<F_m->n_stats; j++){
//      F_stats[j] -= meanstats[j];
        F_stats[j+F_m->n_stats] = F_stats[j];
      }
      F_stats += F_m->n_stats;
//      if (fVerbose){ Rprintf("step %d from %d:\n",i, samplesize);}
      /* This then adds the change statistics to these values */
      if (fVerbose){
        if( ((3*i) % samplesize)==0 && samplesize > 500){
	  Rprintf("Sampled %d from Metropolis-Hastings\n", i);}
      }
      
//      Rprintf("Sampled %d from %d\n", i, samplesize);

    /*********************
    Below is an extremely crude device for letting the user know
    when the chain doesn't accept many of the proposed steps.
    *********************/
//    if (fVerbose){
//      Rprintf("Metropolis-Hastings accepted %7.3f%% of %d steps.\n",
//	      tottaken*100.0/(1.0*interval*samplesize), interval*samplesize); 
//    }
//  }else{
//    if (fVerbose){
//      Rprintf("Metropolis-Hastings accepted %7.3f%% of %d steps.\n",
//	      staken*100.0/(1.0*burnin), burnin); 
//    }
  }
//  Rprintf("netstats: %d\n", samplesize);
//  for (i=0; i < samplesize; i++){
//   for (j=0; j < m->n_stats; j++){
//      Rprintf("%f ", F_stats[j+(m->n_stats)*(i)]);
//   }
//  Rprintf("\n");
//  }
    free(ubar);
    free(u2bar);
}
