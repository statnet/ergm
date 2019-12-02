/*  File src/wtMCMC.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2019 Statnet Commons
 */
#include "wtMCMC.h"
#include "ergm_util.h"
#include "ergm_wtstate.h"
/*****************
 Note on undirected networks:  For j<k, edge {j,k} should be stored
 as (j,k) rather than (k,j).  In other words, only directed networks
 should have (k,j) with k>j.
*****************/

/*****************
 void WtMCMC_wrapper

 Wrapper for a call from R.

 and don't forget that tail -> head
*****************/
SEXP WtMCMC_wrapper(// Network settings
                  SEXP dn, SEXP dflag, SEXP bipartite,
                  // Model settings
                  SEXP nterms, SEXP funnames,
                  SEXP sonames,
                  // Proposal settings
                  SEXP MHProposaltype, SEXP MHProposalpackage,
                  // Numeric inputs
                  SEXP inputs,
                  // Network state
                  SEXP nedges,
                  SEXP tails, SEXP heads,
                  SEXP weights,
                  // MCMC settings
                  SEXP eta, SEXP samplesize, 
                  SEXP burnin, SEXP interval,  
                  SEXP maxedges,
                  SEXP verbose){  
  GetRNGstate();  /* R function enabling uniform RNG */
  ErgmWtState *s = ErgmWtStateInit(// Network settings
                                   asInteger(dn), asInteger(dflag), asInteger(bipartite),
                                   // Model settings
                                   asInteger(nterms), FIRSTCHAR(funnames), FIRSTCHAR(sonames),
                                   // Proposal settings
                                   FIRSTCHAR(MHProposaltype), FIRSTCHAR(MHProposalpackage),
                                   // Numeric inputs
                                   REAL(inputs),
                                   // Network state
                                   asInteger(nedges), (Vertex*) INTEGER(tails), (Vertex*) INTEGER(heads), REAL(weights),
                                   NO_LASTTOGGLE);

  WtNetwork *nwp = s->nwp;
  WtModel *m = s->m;
  WtMHProposal *MHp = s->MHp;

  SEXP sample = PROTECT(allocVector(REALSXP, asInteger(samplesize)*m->n_stats));
  memset(REAL(sample), 0, asInteger(samplesize)*m->n_stats*sizeof(double));

  SEXP status;
  if(MHp) status = PROTECT(ScalarInteger(WtMCMCSample(s,
                                                      REAL(eta), REAL(sample), asInteger(samplesize),
                                                      asInteger(burnin), asInteger(interval), asInteger(maxedges),
                                                      asInteger(verbose))));
  else status = PROTECT(ScalarInteger(MCMC_MH_FAILED));

  const char *outn[] = {"status", "s", WTNWSTATE_NAMES, ""};
  SEXP outl = PROTECT(mkNamed(VECSXP, outn));
  SET_VECTOR_ELT(outl, 0, status);
  SET_VECTOR_ELT(outl, 1, sample);
  
  /* record new generated network to pass back to R */
  if(asInteger(status) == MCMC_OK && asInteger(maxedges)>0){
    WTNWSTATE_SAVE_INTO_RLIST(nwp, outl, 2);
  }

  ErgmWtStateDestroy(s);  
  PutRNGstate();  /* Disable RNG before returning */
  UNPROTECT(3);
  return outl;
}


/*********************
 void WtMCMCSample

 Using the parameters contained in the array eta, obtain the
 network statistics for a sample of size samplesize.  burnin is the
 initial number of Markov chain steps before sampling anything
 and interval is the number of MC steps between successive 
 networks in the sample.  Put all the sampled statistics into
 the networkstatistics array. 
*********************/
MCMCStatus WtMCMCSample(ErgmWtState *s,
			  double *eta, double *networkstatistics, 
			  int samplesize, int burnin, 
			  int interval, int nmax, int verbose) {
  WtNetwork *nwp = s->nwp;
  WtModel *m = s->m;

  int staken, tottaken;
  int i;
    
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
  if(WtMetropolisHastings(s, eta, networkstatistics, burnin, &staken,
			verbose)!=MCMC_OK)
    return MCMC_MH_FAILED;
  if(nmax!=0 && EDGECOUNT(nwp) >= nmax-1){
    ErgmWtStateDestroy(s);  
    error("Number of edges %u exceeds the upper limit set by the user (%u). This can be a sign of degeneracy, but if not, it can be controlled via MCMC.max.maxedges= and/or MCMLE.density.guard= control parameters.", EDGECOUNT(nwp), nmax);
  }
  
/*   if (verbose){ 
       Rprintf(".");
     } */
  
  if (samplesize>1){
    staken = 0;
    tottaken = 0;
    
    /* Now sample networks */
    for (i=1; i < samplesize; i++){
      /* Set current vector of stats equal to previous vector */
      memcpy(networkstatistics+m->n_stats, networkstatistics, m->n_stats*sizeof(double));
      networkstatistics += m->n_stats;
      /* This then adds the change statistics to these values */
      
      /* Catch massive number of edges caused by degeneracy */
      if(WtMetropolisHastings(s, eta, networkstatistics, interval, &staken,
			    verbose)!=MCMC_OK)
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
    if (verbose){
      Rprintf("Sampler accepted %7.3f%% of %lld proposed steps.\n",
	    tottaken*100.0/(1.0*interval*samplesize), (long long) interval*samplesize); 
    }
  }else{
    if (verbose){
      Rprintf("Sampler accepted %7.3f%% of %d proposed steps.\n",
      staken*100.0/(1.0*burnin), burnin); 
    }
  }
  return MCMC_OK;
}

/*********************
 void MetropolisHastings

 In this function, eta is a m->n_stats-vector just as in WtMCMCSample,
 but now networkstatistics is merely another m->n_stats-vector because
 this function merely iterates nsteps times through the Markov
 chain, keeping track of the cumulative change statistics along
 the way, then returns, leaving the updated change statistics in
 the networkstatistics vector.  In other words, this function 
 essentially generates a sample of size one
*********************/
MCMCStatus WtMetropolisHastings (ErgmWtState *s,
				 double *eta, double *networkstatistics,
				 int nsteps, int *staken,
				 int verbose) {
  
  WtNetwork *nwp = s->nwp;
  WtModel *m = s->m;
  WtMHProposal *MHp = s->MHp;

  unsigned int taken=0, unsuccessful=0;
/*  if (verbose)
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
    
    if(verbose>=5){
      Rprintf("MHProposal: ");
      for(unsigned int i=0; i<MHp->ntoggles; i++)
	Rprintf("  (%d, %d) -> %f  ", MHp->toggletail[i], MHp->togglehead[i], MHp->toggleweight[i]);
      Rprintf("\n");
    }

    /* Calculate change statistics,
       remembering that tail -> head */
    WtChangeStats(MHp->ntoggles, MHp->toggletail, MHp->togglehead, MHp->toggleweight, nwp, m);

    if(verbose>=5){
      Rprintf("Changes: (");
      for(unsigned int i=0; i<m->n_stats; i++)
	Rprintf(" %f ", m->workspace[i]);
      Rprintf(")\n");
    }
    
    /* Calculate inner (dot) product */
    double ip = dotprod(eta, m->workspace, m->n_stats);

    /* The logic is to set cutoff = ip+logratio ,
       then let the MH probability equal min{exp(cutoff), 1.0}.
       But we'll do it in log space instead.  */
    double cutoff = ip + MHp->logratio;

    if(verbose>=5){
      Rprintf("log acceptance probability: %f + %f = %f\n", ip, MHp->logratio, cutoff);
    }
    
    /* if we accept the proposed network */
    if (cutoff >= 0.0 || logf(unif_rand()) < cutoff) { 
      if(verbose>=5){
	Rprintf("Accepted.\n");
      }

      /* Make proposed toggles (updating timestamps--i.e., for real this time) */
      for(unsigned int i=0; i < MHp->ntoggles; i++){
	Vertex t=MHp->toggletail[i], h=MHp->togglehead[i];
	double w=MHp->toggleweight[i];

	WtGET_EDGE_UPDATE_STORAGE_SET(t, h, w, nwp, m, MHp);
      }
      /* record network statistics for posterity */
      for (unsigned int i = 0; i < m->n_stats; i++){
	networkstatistics[i] += m->workspace[i];
      }
      taken++;
    }else{
      if(verbose>=5){
	Rprintf("Rejected.\n");
      }
    }
  }
  
  *staken = taken;
  return MCMC_OK;
}

