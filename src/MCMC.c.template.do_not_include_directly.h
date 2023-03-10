/*  File src/MCMC.c.template.do_not_include_directly.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2023 Statnet Commons
 */
#include "ergm_etamap.h"
#include "ergm_util.h"
/*****************
 Note on undirected networks:  For j<k, edge {j,k} should be stored
 as (j,k) rather than (k,j).  In other words, only directed networks
 should have (k,j) with k>j.
*****************/

/*****************
 void DISPATCH_MCMC_wrapper

 Wrapper for a call from R.

 and don't forget that tail -> head
*****************/
SEXP DISPATCH_MCMC_wrapper(SEXP stateR,
                    // MCMC settings
                    SEXP eta, SEXP samplesize,SEXP burnin, SEXP interval,SEXP maxedges,SEXP verbose){
  GetRNGstate();  /* R function enabling uniform RNG */
  unsigned int protected = 0;

  DISPATCH_ErgmState *s = DISPATCH_ErgmStateInit(stateR, 0);
  if(asInteger(maxedges) < 0){
    s->save = PROTECT(allocVector(VECSXP, asInteger(samplesize))); protected++;
  }else s->save = NULL;

  DISPATCH_Model *m = s->m;
  DISPATCH_MHProposal *MHp = s->MHp;

  SEXP sample = PROTECT(allocVector(REALSXP, asInteger(samplesize)*m->n_stats)); protected++;
  memset(REAL(sample), 0, asInteger(samplesize)*m->n_stats*sizeof(double));
  memcpy(REAL(sample), s->stats, m->n_stats*sizeof(double));

  SEXP status;
  if(MHp) status = PROTECT(ScalarInteger(DISPATCH_MCMCSample(s,
                                                             REAL(eta), REAL(sample), asInteger(samplesize),
                                                             asInteger(burnin), asInteger(interval), abs(asInteger(maxedges)),
                                                             asInteger(verbose))));
  else status = PROTECT(ScalarInteger(MCMC_MH_FAILED));
  protected++;

  const char *outn[] = {"status", "s", "state", "saved", ""};
  SEXP outl = PROTECT(mkNamed(VECSXP, outn)); protected++;
  SET_VECTOR_ELT(outl, 0, status);
  SET_VECTOR_ELT(outl, 1, sample);

  /* record new generated network to pass back to R */
  if(asInteger(status) == MCMC_OK && asInteger(maxedges) != 0){
    s->stats = REAL(sample) + (asInteger(samplesize)-1)*m->n_stats;
    SET_VECTOR_ELT(outl, 2, DISPATCH_ErgmStateRSave(s));
  }

  if(s->save) SET_VECTOR_ELT(outl, 3, s->save);

  DISPATCH_ErgmStateDestroy(s);
  PutRNGstate();  /* Disable RNG before returning */
  UNPROTECT(protected); protected = 0;
  return outl;
}


/*********************
 void DISPATCH_MCMCSample

 Using the parameters contained in the array eta, obtain the
 network statistics for a sample of size samplesize.  burnin is the
 initial number of Markov chain steps before sampling anything
 and interval is the number of MC steps between successive
 networks in the sample.  Put all the sampled statistics into
 the networkstatistics array.
*********************/
MCMCStatus DISPATCH_MCMCSample(DISPATCH_ErgmState *s,
			  double *eta, double *networkstatistics,
			  int samplesize, int burnin,
			  int interval, int nmax, int verbose) {
  DISPATCH_Network *nwp = s->nwp;
  DISPATCH_Model *m = s->m;

  int staken, tottaken;

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
  if(DISPATCH_MetropolisHastings(s, eta, networkstatistics, burnin, &staken,verbose)!=MCMC_OK)
    return MCMC_MH_FAILED;
  if(nmax!=0 && EDGECOUNT(nwp) >= nmax-1){
    return MCMC_TOO_MANY_EDGES;
  }

  if(s->save){
    s->stats = networkstatistics;
    SET_VECTOR_ELT(s->save, 0, DISPATCH_ErgmStateRSave(s));
  }

/*   if (verbose){
       Rprintf(".");
     } */

  if (samplesize>1){
    staken = 0;
    tottaken = 0;

    /* Now sample networks */
    for (unsigned int i=1; i < samplesize; i++){
      /* Set current vector of stats equal to previous vector */
      memcpy(networkstatistics+m->n_stats, networkstatistics, m->n_stats*sizeof(double));
      networkstatistics += m->n_stats;
      /* This then adds the change statistics to these values */

      /* Catch massive number of edges caused by degeneracy */
      if(DISPATCH_MetropolisHastings(s, eta, networkstatistics, interval, &staken,
			    verbose)!=MCMC_OK)
	return MCMC_MH_FAILED;
      if(nmax!=0 && EDGECOUNT(nwp) >= nmax-1){
	return MCMC_TOO_MANY_EDGES;
      }

      if(s->save){
        s->stats = networkstatistics;
        SET_VECTOR_ELT(s->save, i, DISPATCH_ErgmStateRSave(s));
      }

      tottaken += staken;

      R_CheckUserInterrupt();
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

 In this function, eta is a m->n_stats-vector just as in DISPATCH_MCMCSample,
 but now networkstatistics is merely another m->n_stats-vector because
 this function merely iterates nsteps times through the Markov
 chain, keeping track of the cumulative change statistics along
 the way, then returns, leaving the updated change statistics in
 the networkstatistics vector.  In other words, this function
 essentially generates a sample of size one
*********************/
MCMCStatus DISPATCH_MetropolisHastings (DISPATCH_ErgmState *s,
				 double *eta, double *networkstatistics,
				 int nsteps, int *staken,
				 int verbose) {

  DISPATCH_Network *nwp = s->nwp;
  DISPATCH_Model *m = s->m;
  DISPATCH_MHProposal *MHp = s->MHp;

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
      for(unsigned int i=0; i<MHp->ntoggles; i++) PROP_PRINT;
      Rprintf("\n");
    }

    /* Calculate change statistics,
       remembering that tail -> head */
    PROP_CHANGESTATS;

    if(verbose>=5) print_vector("stat diff", m->workspace, m->n_stats);

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
      for(unsigned int i=0; i < MHp->ntoggles; i++) PROP_COMMIT;
      /* record network statistics for posterity */
      addonto(networkstatistics, m->workspace, m->n_stats);
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

/* *** don't forget tail -> head */

SEXP DISPATCH_MCMCPhase12 (SEXP stateR,
                           // Phase12 settings
                           SEXP theta0,
                           SEXP burnin, SEXP interval,
                           SEXP gain, SEXP phase1, SEXP nsub,
                           SEXP min_iterations, SEXP max_iterations,
                           SEXP maxedges,
                           SEXP verbose){
  GetRNGstate();  /* R function enabling uniform RNG */
  DISPATCH_ErgmState *s = DISPATCH_ErgmStateInit(stateR, 0);

  DISPATCH_Model *m = s->m;
  DISPATCH_MHProposal *MHp = s->MHp;

  SEXP stats = PROTECT(allocVector(REALSXP, m->n_stats));
  memcpy(REAL(stats), s->stats, m->n_stats*sizeof(double));
  s->stats = REAL(stats);

  unsigned int n_param = length(theta0);
  SEXP theta = PROTECT(allocVector(REALSXP, n_param));
  memcpy(REAL(theta), REAL(theta0), n_param*sizeof(double));

  SEXP status;
  if(MHp) status = PROTECT(ScalarInteger(DISPATCH_MCMCSamplePhase12(s,
                                                           REAL(theta), n_param, asReal(gain), asInteger(phase1), asInteger(nsub),
                                                           asInteger(min_iterations), asInteger(max_iterations),
                                                           asInteger(burnin), asInteger(interval),
                                                           asInteger(verbose))));
  else status = PROTECT(ScalarInteger(MCMC_MH_FAILED));

  const char *outn[] = {"status", "theta", "state", ""};
  SEXP outl = PROTECT(mkNamed(VECSXP, outn));
  SET_VECTOR_ELT(outl, 0, status);
  SET_VECTOR_ELT(outl, 1, theta);

  /* record new generated network to pass back to R */
  if(asInteger(status) == MCMC_OK && asInteger(maxedges)>0){
    SET_VECTOR_ELT(outl, 2, DISPATCH_ErgmStateRSave(s));
  }

  DISPATCH_ErgmStateDestroy(s);
  PutRNGstate();  /* Disable RNG before returning */
  UNPROTECT(4);
  return outl;
}

/*********************
 void MCMCSamplePhase12
*********************/
MCMCStatus DISPATCH_MCMCSamplePhase12(DISPATCH_ErgmState *s,
                               double *theta, unsigned int n_param, double gain, int nphase1, int nsubphases,
                               int min_iterations, int max_iterations,
                               int burnin,
                               int interval, int verbose){
  DISPATCH_Model *m = s->m;
  SEXP etamap = getListElement(m->R, "etamap");
  const int *theta_offset = LOGICAL(getListElement(etamap, "offsettheta"));
  const double *theta_min = REAL(getListElement(etamap, "mintheta"));
  const double *theta_max = REAL(getListElement(etamap, "maxtheta"));
  const unsigned int n_stats = m->n_stats;

  int staken;

  /*Rprintf("nsubphases %d\n", nsubphases); */

  double *ubar = R_calloc(n_param, double),
    *u2bar = R_calloc(n_param, double),
    *aDdiaginv = R_calloc(n_param, double);

  /*********************
   Burn in step.  While we're at it, use burnin statistics to
   prepare covariance matrix for Mahalanobis distance calculations
   in subsequent calls to M-H
   *********************/

  double *eta = R_calloc(n_stats, double),
    *etagrad = R_calloc(n_stats*n_param, double);

  ergm_eta(theta, etamap, eta);
  ergm_etagrad(theta, etamap, etagrad);

  staken = 0;
  Rprintf("Starting burnin of %d steps\n", burnin);
  MCMCStatus status = DISPATCH_MetropolisHastings(s, eta,
                                           s->stats, burnin, &staken,
                                           verbose);
  if(status!=MCMC_OK) return status;
  Rprintf("Phase 1: %d steps (interval = %d)\n", nphase1,interval);
  /* Now sample networks */
  for (unsigned int i=0; i <= nphase1; i++){
    MCMCStatus status = DISPATCH_MetropolisHastings(s, eta,
                                                    s->stats, interval, &staken,
                                                    verbose);
    if(status!=MCMC_OK) return status;
    if(i > 0){
      for(unsigned int j=0; j<n_param; j++){
        double u = 0;
        for(unsigned int k=0; k<n_stats; k++) u += etagrad[j+n_param*k] * s->stats[k];
        ubar[j]  += u;
        u2bar[j] += u*u;
      }
    }
  }

  if (verbose){
    Rprintf("Returned from Phase 1\n");
  }
  for(unsigned int j=0; j<n_param; j++){
    if(theta_offset[j]) continue;
    aDdiaginv[j] = u2bar[j] - ubar[j]*ubar[j]/nphase1;
    if(aDdiaginv[j] > 0.0){
      aDdiaginv[j] = nphase1*gain/aDdiaginv[j];
    }else{
      aDdiaginv[j]=0.00001;
    }
  }

  if(verbose) print_vector("gain * V^-1", aDdiaginv, n_param);

  if(verbose){
    Rprintf("\nPhase 2:\n");
  }

  double *theta_sum = R_calloc(n_param, double),
    *esteq = R_calloc(n_param, double),
    *esteq_old = R_calloc(n_param, double),
    *esteq_prod_cum = R_calloc(n_param, double);

  /* Now run Phase2, In which there are several subphases. We sample networks in each subphases */
  for(unsigned int subphase = 1; subphase <= nsubphases; subphase++){
    unsigned int N2klower = trunc(min_iterations*pow(2.52,(subphase-1)))+1; /* The lower bound for the number of iterations in subphase. */
    unsigned int N2kupper = max_iterations - min_iterations + N2klower; /* The Upper bound for the number of iterations in subphase. */

    memset(theta_sum, 0, n_param*sizeof(double));
    memset(esteq_prod_cum, 0, n_param*sizeof(double));

    for(unsigned int i=1; ; i++){
      MCMCStatus status = DISPATCH_MetropolisHastings(s, eta, s->stats, interval, &staken,verbose); /*Take a sample network*/
      if(status!=MCMC_OK) return status;

      memcpy(esteq_old, esteq, n_param*sizeof(double));
      memset(esteq, 0, n_param*sizeof(double));

      if(verbose>=4){
        Rprintf("\n");
        print_vector("eta", eta, n_stats);
        print_vector("statistic", s->stats, n_stats);
        print_matrix("Dtheta/Deta", etagrad, n_param, n_stats);
      }

      for (unsigned int j=0; j<n_param; j++){
        if(theta_offset[j]) continue;

        for(unsigned int k=0; k<n_stats; k++)
          esteq[j] += etagrad[j+n_param*k] * s->stats[k];

        esteq_prod_cum[j] += esteq[j] * esteq_old[j];

        /*Update Theta*/
        theta[j] -= aDdiaginv[j] * esteq[j];
        if(theta[j] < theta_min[j]) theta[j] = theta_min[j];
        if(theta[j] > theta_max[j]) theta[j] = theta_max[j];

        /*Take summation of theta for average calculation*/
        theta_sum[j] += theta[j];
      }

      ergm_eta(theta, etamap, eta);
      ergm_etagrad(theta, etamap, etagrad);
      if(verbose>=4){
        print_vector("estimating function", esteq, n_param);
        print_vector("successive change products", esteq_prod_cum, n_param);
        print_vector("theta", theta, n_param);
      }

      Rboolean subphase_done = FALSE;

      if(i >= N2kupper){
        subphase_done = TRUE;
        if(verbose>=3) Rprintf("Subphase ran out of steps.\n");
      }else if(i >= N2klower){
        subphase_done = TRUE;
        for(unsigned int j=0; j<n_param; j++){
          if(!theta_offset[j] && esteq_prod_cum[j] >= 0){
            subphase_done = FALSE;
            break;
          }
        }
        if(subphase_done && verbose>=3) Rprintf("Suphase reached the stopping criteria.\n");
      }

      if(subphase_done){
        for(unsigned int j=0; j<n_param; j++){
          if(theta_offset[j]) continue;
          if(subphase == nsubphases) theta[j] = theta_sum[j] / i;
          aDdiaginv[j] /= 2.0;
        }
        break;
      }
    }

    /* Set current vector of stats equal to previous vector */
    R_CheckUserInterrupt();
#ifdef Win32
    R_FlushConsole();
    R_ProcessEvents();
#endif
  }

  // R_calloc()-ed variables freed on return to R.

  return MCMC_OK;
}
