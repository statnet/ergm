/*  File src/MCMC.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2019 Statnet Commons
 */
#include "MCMC.h"
#include "ergm_util.h"
/*****************
 Note on undirected networks:  For j<k, edge {j,k} should be stored
 as (j,k) rather than (k,j).  In other words, only directed networks
 should have (k,j) with k>j.
*****************/

#define TOINTSXP(x) x = PROTECT(coerceVector(x, INTSXP))
#define TOREALSXP(x) x = PROTECT(coerceVector(x, REALSXP))
#define FIRSTCHAR(x) CHAR(STRING_ELT(x, 0))

/*****************
 void MCMC_wrapper

 Wrapper for a call from R.

 and don't forget that tail -> head
*****************/
SEXP MCMC_wrapper(// Network settings
                  SEXP dn, SEXP dflag, SEXP bipartite,
                  // Model settings
                  SEXP nterms, SEXP funnames,
                  SEXP sonames,
                  // Proposal settings
                  SEXP MHProposaltype, SEXP MHProposalpackage,
                  SEXP attribs, SEXP maxout, SEXP maxin, SEXP minout,
                  SEXP minin, SEXP condAllDegExact, SEXP attriblength,
                  // Numeric inputs
                  SEXP inputs,
                  // Network state
                  SEXP nedges,
                  SEXP tails, SEXP heads,
                  // MCMC settings
                  SEXP eta, SEXP samplesize, 
                  SEXP burnin, SEXP interval,  
                  SEXP maxedges,
                  SEXP verbose){  
  GetRNGstate();  /* R function enabling uniform RNG */
  // Ensure correct data types
  TOINTSXP(tails);
  TOINTSXP(heads);
  TOREALSXP(inputs);
  TOINTSXP(attribs);
  TOINTSXP(maxout);
  TOINTSXP(maxin);
  TOINTSXP(minout);
  TOINTSXP(minin);
  ErgmState *s = ErgmStateInit(// Network settings
                               asInteger(dn), asInteger(dflag), asInteger(bipartite),
                               // Model settings
                               asInteger(nterms), FIRSTCHAR(funnames), FIRSTCHAR(sonames),
                               // Proposal settings
                               FIRSTCHAR(MHProposaltype), FIRSTCHAR(MHProposalpackage), INTEGER(attribs), INTEGER(maxout), INTEGER(maxin), INTEGER(minout), INTEGER(minin), asInteger(condAllDegExact), asInteger(attriblength),
                               // Numeric inputs
                               REAL(inputs),
                               // Network state
                               asInteger(nedges), (Vertex*) INTEGER(tails), (Vertex*) INTEGER(heads));
  UNPROTECT(9);

  Network *nwp = s->nwp;
  Model *m = s->m;
  MHProposal *MHp = s->MHp;

  SEXP sample = PROTECT(allocVector(REALSXP, asInteger(samplesize)*m->n_stats));
  
  TOREALSXP(eta);
  SEXP status;
  if(MHp) status = PROTECT(ScalarInteger(MCMCSample(s,
                                                    REAL(eta), REAL(sample), asInteger(samplesize),
                                                    asInteger(burnin), asInteger(interval),
                                                    asInteger(verbose), asInteger(maxedges))));
  else status = PROTECT(ScalarInteger(MCMC_MH_FAILED));

  SEXP outl = PROTECT(allocVector(VECSXP, 4));
  SET_VECTOR_ELT(outl, 0, status);
  SET_VECTOR_ELT(outl, 1, sample);
  
  /* record new generated network to pass back to R */
  if(asInteger(status) == MCMC_OK && asInteger(maxedges)>0){
    SEXP newnetworktails, newnetworkheads;
    newnetworktails = PROTECT(allocVector(INTSXP, EDGECOUNT(nwp)+1));
    newnetworkheads = PROTECT(allocVector(INTSXP, EDGECOUNT(nwp)+1));

    INTEGER(newnetworktails)[0]=INTEGER(newnetworkheads)[0]=
      EdgeTree2EdgeList((Vertex*)INTEGER(newnetworktails)+1,
			(Vertex*)INTEGER(newnetworkheads)+1,
			nwp,asInteger(maxedges)-1);

    SET_VECTOR_ELT(outl, 2, newnetworktails);
    SET_VECTOR_ELT(outl, 3, newnetworkheads);
    UNPROTECT(2);
  }

  ErgmStateDestroy(s);  
  PutRNGstate();  /* Disable RNG before returning */
  UNPROTECT(5);
  return outl;
}


/*********************
 MCMCStatus MCMCSample

 Using the parameters contained in the array theta, obtain the
 network statistics for a sample of size samplesize.  burnin is the
 initial number of Markov chain steps before sampling anything
 and interval is the number of MC steps between successive 
 networks in the sample.  Put all the sampled statistics into
 the networkstatistics array. 
*********************/
MCMCStatus MCMCSample(ErgmState *s,
                      double *theta, double *networkstatistics, 
                      int samplesize, int burnin, 
                      int interval, int fVerbose, int nmax){
  Network *nwp = s->nwp;
  Model *m = s->m;
  // MHProposal *MHp = s->MHp;

  int staken, tottaken;
  int i, j;
    
  /*********************
  networkstatistics are modified in groups of m->n_stats, and they
  reflect the CHANGE in the values of the statistics from the
  original (observed) network.  Thus, when we begin, the initial 
  values of the first group of m->n_stats networkstatistics should 
  all be zero
  *********************/
/*for (j=0; j < m->n_stats; j++) */
/*  networkstatistics[j] = 0.0; */
/* Rprintf("\n"); */
/* for (j=0; j < m->n_stats; j++){ */
/*   Rprintf("j %d %f\n",j,networkstatistics[j]); */
/* } */
/* Rprintf("\n"); */

  /*********************
   Burn in step.
   *********************/
/*  Catch more edges than we can return */
  if(MetropolisHastings(s, theta, networkstatistics, burnin, &staken,
			fVerbose)!=MCMC_OK)
    return MCMC_MH_FAILED;
  if(nmax!=0 && EDGECOUNT(nwp) >= nmax-1){
    ErgmStateDestroy(s);  
    error("Number of edges %u exceeds the upper limit set by the user (%u). This can be a sign of degeneracy, but if not, it can be controlled via MCMC.max.maxedges= and/or MCMLE.density.guard= control parameters.", EDGECOUNT(nwp), nmax);
  }
  
/*   if (fVerbose){ 
       Rprintf(".");
     } */
  
  if (samplesize>1){
    staken = 0;
    tottaken = 0;
    
    /* Now sample networks */
    for (i=1; i < samplesize; i++){
      /* Set current vector of stats equal to previous vector */
      for (j=0; j<m->n_stats; j++){
        networkstatistics[j+m->n_stats] = networkstatistics[j];
      }
      networkstatistics += m->n_stats;
      /* This then adds the change statistics to these values */
      
      /* Catch massive number of edges caused by degeneracy */
      if(MetropolisHastings(s, theta, networkstatistics, interval, &staken,
			    fVerbose)!=MCMC_OK)
	return MCMC_MH_FAILED;
      if(nmax!=0 && EDGECOUNT(nwp) >= nmax-1){
	return MCMC_TOO_MANY_EDGES;
      }
      tottaken += staken;

#ifdef Win32
      if( ((100*i) % samplesize)==0 && samplesize > 500){
	R_FlushConsole();
    	R_ProcessEvents();
      }
#endif
    }
    /*********************
    Below is an extremely crude device for letting the user know
    when the chain doesn't accept many of the proposed steps.
    *********************/
    if (fVerbose){
      Rprintf("Sampler accepted %7.3f%% of %lld proposed steps.\n",
	    tottaken*100.0/(1.0*interval*samplesize), (long long) interval*samplesize); 
    }
  }else{
    if (fVerbose){
      Rprintf("Sampler accepted %7.3f%% of %d proposed steps.\n",
      staken*100.0/(1.0*burnin), burnin); 
    }
  }
  return MCMC_OK;
}

/*********************
 void MetropolisHastings

 In this function, theta is a m->n_stats-vector just as in MCMCSample,
 but now networkstatistics is merely another m->n_stats-vector because
 this function merely iterates nsteps times through the Markov
 chain, keeping track of the cumulative change statistics along
 the way, then returns, leaving the updated change statistics in
 the networkstatistics vector.  In other words, this function 
 essentially generates a sample of size one
*********************/
MCMCStatus MetropolisHastings(ErgmState *s,
			      double *theta, double *networkstatistics,
			      int nsteps, int *staken,
			      int fVerbose) {

  Network *nwp = s->nwp;
  Model *m = s->m;
  MHProposal *MHp = s->MHp;

  unsigned int taken=0, unsuccessful=0;
/*  if (fVerbose)
    Rprintf("Now proposing %d MH steps... ", nsteps); */
  for(unsigned int step=0; step < nsteps; step++) {
    MHp->logratio = 0;
    (*(MHp->p_func))(MHp, nwp); /* Call MH function to propose toggles */

    if(MHp->toggletail[0]==MH_FAILED){
      switch(MHp->togglehead[0]){
      case MH_UNRECOVERABLE:
	error("Something very bad happened during proposal. Memory has not been deallocated, so restart R soon.");
	
      case MH_IMPOSSIBLE:
	Rprintf("MH MHProposal function encountered a configuration from which no toggle(s) can be proposed.\n");
	return MCMC_MH_FAILED;
	
      case MH_UNSUCCESSFUL:
	warning("MH MHProposal function failed to find a valid proposal.");
	unsuccessful++;
	if(unsuccessful>taken*MH_QUIT_UNSUCCESSFUL){
	  Rprintf("Too many MH MHProposal function failures.\n");
	  return MCMC_MH_FAILED;
	}
      case MH_CONSTRAINT:
	continue;
      }
    }
    
    if(fVerbose>=5){
      Rprintf("MHProposal: ");
      for(unsigned int i=0; i<MHp->ntoggles; i++)
	Rprintf(" (%d, %d)", MHp->toggletail[i], MHp->togglehead[i]);
      Rprintf("\n");
    }

    /* Calculate change statistics,
       remembering that tail -> head */
    ChangeStats(MHp->ntoggles, MHp->toggletail, MHp->togglehead, nwp, m);

    if(fVerbose>=5){
      Rprintf("Changes: (");
      for(unsigned int i=0; i<m->n_stats; i++)
	Rprintf(" %f ", m->workspace[i]);
      Rprintf(")\n");
    }
    
    /* Calculate inner (dot) product */
    double ip = dotprod(theta, m->workspace, m->n_stats);

    /* The logic is to set cutoff = ip+logratio ,
       then let the MH probability equal min{exp(cutoff), 1.0}.
       But we'll do it in log space instead.  */
    double cutoff = ip + MHp->logratio;

    if(fVerbose>=5){
      Rprintf("log acceptance probability: %f + %f = %f\n", ip, MHp->logratio, cutoff);
    }
    
    /* if we accept the proposed network */
    if (cutoff >= 0.0 || logf(unif_rand()) < cutoff) { 
      if(fVerbose>=5){
	Rprintf("Accepted.\n");
      }

      /* Make proposed toggles (updating timestamps--i.e., for real this time) */
      for(unsigned int i=0; i < MHp->ntoggles; i++){
	GET_EDGE_UPDATE_STORAGE_TOGGLE(MHp->toggletail[i], MHp->togglehead[i], nwp, m, MHp);
      }
      /* record network statistics for posterity */
      for (unsigned int i = 0; i < m->n_stats; i++){
	networkstatistics[i] += m->workspace[i];
      }
      taken++;
    }else{
      if(fVerbose>=5){
	Rprintf("Rejected.\n");
      }
    }
  }
  
  *staken = taken;
  return MCMC_OK;
}

/* *** don't forget tail -> head */

void MCMCPhase12 (int *tails, int *heads, int *dnedges, 
		  int *dn, int *dflag, int *bipartite, 
		  int *nterms, char **funnames,
		  char **sonames, 
		  char **MHProposaltype, char **MHProposalpackage,
		  double *inputs, 
		  double *theta0, int *samplesize,
		  double *gain, double *meanstats, int *phase1, int *nsub,
		  double *sample, int *burnin, int *interval,  
		  int *newnetworktails, 
		  int *newnetworkheads, 
		  int *fVerbose, 
		  int *attribs, int *maxout, int *maxin, int *minout,
		  int *minin, int *condAllDegExact, int *attriblength, 
		  int *maxedges,
		  int *mtails, int *mheads, int *mdnedges)  {
  int directed_flag;
  int nphase1, nsubphases;
  Vertex n_nodes, bip;
  Edge n_edges, nmax;
  
  nphase1 = *phase1; 
  nsubphases = *nsub;

  n_nodes = (Vertex)*dn; 
  n_edges = (Edge)*dnedges; 
  nmax = (Edge)*maxedges; 
  bip = (Vertex)*bipartite; 
  
  GetRNGstate();  /* R function enabling uniform RNG */
  
  directed_flag = *dflag;

  ErgmState *s = ErgmStateInit(// Network settings
                               n_nodes, directed_flag, bip,
                               // Model settings
                               *nterms, *funnames, *sonames,
                               // Proposal settings
                               *MHProposaltype, *MHProposalpackage, attribs, maxout, maxin, minout, minin, *condAllDegExact, *attriblength,
                               // Numeric inputs
                               inputs,
                               // Network state
                               n_edges, (Vertex*) tails, (Vertex*) heads);

  Network *nwp = s->nwp;
  MHProposal *MHp = s->MHp;

  if(MHp)
    MCMCSamplePhase12(s,
		      theta0, *gain, meanstats, nphase1, nsubphases, sample, *samplesize,
		      *burnin, *interval,
		      (int)*fVerbose);
  else error("MH Proposal failed.");

  /* record new generated network to pass back to R */
  if(nmax>0 && newnetworktails && newnetworkheads)
    newnetworktails[0]=newnetworkheads[0]=EdgeTree2EdgeList((Vertex*)newnetworktails+1,(Vertex*)newnetworkheads+1,nwp,nmax-1);

  ErgmStateDestroy(s);
  PutRNGstate();  /* Disable RNG before returning */
}

/*********************
 void MCMCSamplePhase12

 Using the parameters contained in the array theta, obtain the
 network statistics for a sample of size samplesize.  burnin is the
 initial number of Markov chain steps before sampling anything
 and interval is the number of MC steps between successive 
 networks in the sample.  Put all the sampled statistics into
 the networkstatistics array. 
*********************/
void MCMCSamplePhase12(ErgmState *s,
		       double *theta, double gain, double *meanstats, int nphase1, int nsubphases, double *networkstatistics, 
		       int samplesize, int burnin, 
		       int interval, int fVerbose){
  Model *m = s->m;

  int staken, tottaken, ptottaken;
  int i, j, iter=0;
  
/*Rprintf("nsubphases %d\n", nsubphases); */

  /*if (fVerbose)
    Rprintf("The number of statistics is %i and the total samplesize is %d\n",
    m->n_stats,samplesize);*/

  /*********************
  networkstatistics are modified in groups of m->n_stats, and they
  reflect the CHANGE in the values of the statistics from the
  original (observed) network.  Thus, when we begin, the initial 
  values of the first group of m->n_stats networkstatistics should 
  all be zero
  *********************/
  double *ubar, *u2bar, *aDdiaginv;
  ubar = (double *)Calloc(m->n_stats, double);
  u2bar = (double *)Calloc(m->n_stats, double);
  aDdiaginv = (double *)Calloc(m->n_stats, double);
  for (j=0; j < m->n_stats; j++){
    networkstatistics[j] = -meanstats[j];
    ubar[j] = 0.0;
    u2bar[j] = 0.0;
  }

  /*********************
   Burn in step.  While we're at it, use burnin statistics to 
   prepare covariance matrix for Mahalanobis distance calculations 
   in subsequent calls to M-H
   *********************/
/*Rprintf("MCMCSampleDyn pre burnin numdissolve %d\n", *numdissolve); */
  
    staken = 0;
    Rprintf("Starting burnin of %d steps\n", burnin);
    MetropolisHastings(s, theta,
                       networkstatistics, burnin, &staken,
                       fVerbose);
    Rprintf("Phase 1: %d steps (interval = %d)\n", nphase1,interval);
    /* Now sample networks */
    for (i=0; i <= nphase1; i++){
      MetropolisHastings(s, theta,
                         networkstatistics, interval, &staken,
                         fVerbose);
      if(i > 0){
       for (j=0; j<m->n_stats; j++){
        ubar[j]  += networkstatistics[j];
        u2bar[j] += networkstatistics[j]*networkstatistics[j];
/*  Rprintf("j %d ubar %f u2bar %f ns %f\n", j,  ubar[j], u2bar[j], */
/*		  networkstatistics[j]); */
       }
      }
    }
    if (fVerbose){
      Rprintf("Returned from Phase 1\n");
      Rprintf("\n gain times inverse variances:\n");
    }
    for (j=0; j<m->n_stats; j++){
      aDdiaginv[j] = u2bar[j]-ubar[j]*ubar[j]/(1.0*nphase1);
      if( aDdiaginv[j] > 0.0){
        aDdiaginv[j] = nphase1*gain/aDdiaginv[j];
      }else{
	aDdiaginv[j]=0.00001;
      }
      if (fVerbose){ Rprintf(" %f", aDdiaginv[j]);}
    }
    if (fVerbose){ Rprintf("\n"); }
  
    staken = 0;
    tottaken = 0;
    ptottaken = 0;
    
    if (fVerbose){
      Rprintf("Phase 2: (samplesize = %d)\n", samplesize);
    }
    /* Now sample networks */
    for (i=1; i < samplesize; i++){
      
      MetropolisHastings(s, theta,
                         networkstatistics, interval, &staken,
                         fVerbose);
    /* Update theta0 */
/*Rprintf("initial:\n"); */
      for (j=0; j<m->n_stats; j++){
        theta[j] -= aDdiaginv[j] * networkstatistics[j];
      }
/*Rprintf("\n"); */
/*    if (fVerbose){ Rprintf("nsubphases %d i %d\n", nsubphases, i); } */
      if (i==(nsubphases)){
	nsubphases = trunc(nsubphases*2.52) + 1;
        if (fVerbose){
	 iter++;
	 Rprintf("End of iteration %d; Updating the number of sub-phases to be %d\n",iter,nsubphases);
	}
        for (j=0; j<m->n_stats; j++){
          aDdiaginv[j] /= 2.0;
          if (fVerbose){Rprintf("theta_%d = %f; change statistic[%d] = %f\n",
		                 j+1, theta[j], j+1, networkstatistics[j]);}
/*        if (fVerbose){ Rprintf(" %f statsmean %f",  theta[j],(networkstatistics[j]-meanstats[j])); } */
        }
        if (fVerbose){ Rprintf("\n"); }
      }
      /* Set current vector of stats equal to previous vector */
      for (j=0; j<m->n_stats; j++){
/*      networkstatistics[j] -= meanstats[j]; */
        networkstatistics[j+m->n_stats] = networkstatistics[j];
      }
      networkstatistics += m->n_stats;
/*      if (fVerbose){ Rprintf("step %d from %d:\n",i, samplesize);} */
      /* This then adds the change statistics to these values */
      tottaken += staken;
#ifdef Win32
      if( ((100*i) % samplesize)==0 && samplesize > 500){
	R_FlushConsole();
    	R_ProcessEvents();
      }
#endif
      if (fVerbose){
        if( ((3*i) % samplesize)==0 && samplesize > 500){
        Rprintf("Sampled %d from Metropolis-Hastings\n", i);}
      }
      
      if( ((3*i) % samplesize)==0 && tottaken == ptottaken){
        ptottaken = tottaken; 
        Rprintf("Warning:  Metropolis-Hastings algorithm has accepted only "
        "%d steps out of a possible %d\n",  ptottaken-tottaken, i); 
      }
/*      Rprintf("Sampled %d from %d\n", i, samplesize); */

    /*********************
    Below is an extremely crude device for letting the user know
    when the chain doesn't accept many of the proposed steps.
    *********************/
/*    if (fVerbose){ */
/*      Rprintf("Metropolis-Hastings accepted %7.3f%% of %d steps.\n", */
/*	      tottaken*100.0/(1.0*interval*samplesize), interval*samplesize);  */
/*    } */
/*  }else{ */
/*    if (fVerbose){ */
/*      Rprintf("Metropolis-Hastings accepted %7.3f%% of %d steps.\n", */
/*	      staken*100.0/(1.0*burnin), burnin);  */
/*    } */
  }
/*  Rprintf("netstats: %d\n", samplesize); */
/*  for (i=0; i < samplesize; i++){ */
/*   for (j=0; j < m->n_stats; j++){ */
/*      Rprintf("%f ", networkstatistics[j+(m->n_stats)*(i)]); */
/*   } */
/*  Rprintf("\n"); */
/*  } */
  if (fVerbose){
    Rprintf("Phase 3: MCMC-Newton-Raphson\n");
  }

  Free(ubar);
  Free(u2bar);
}

