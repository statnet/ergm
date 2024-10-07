/*  File src/SAN.c.template.do_not_include_directly.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2024 Statnet Commons
 */

/*****************
 Note on undirected networks:  For j<k, edge {j,k} should be stored
 as (j,k) rather than (k,j).  In other words, only directed networks
 should have (k,j) with k>j.
*****************/

/*****************
 void DISPATCH_SAN_wrapper

 Wrapper for a call from R.
*****************/

SEXP DISPATCH_SAN_wrapper(SEXP stateR,
                   // MCMC settings
                   SEXP tau,
                   SEXP samplesize, SEXP nsteps,
                   SEXP invcov,
                   SEXP statindices,
                   SEXP offsetindices,
                   SEXP offsets,
                   SEXP verbose){
  GetRNGstate();  /* R function enabling uniform RNG */
  unsigned int nstats = length(statindices), noffsets = length(offsetindices);
  
  DISPATCH_ErgmState *s = DISPATCH_ErgmStateInit(stateR, 0);
  DISPATCH_MHProposal *MHp = s->MHp;

  SEXP sample = PROTECT(allocVector(REALSXP, asInteger(samplesize)*nstats));
  memset(REAL(sample), 0, asInteger(samplesize)*nstats*sizeof(double));
  memcpy(REAL(sample), s->stats, nstats*sizeof(double));
  SEXP prop_sample = PROTECT(allocVector(REALSXP, asInteger(samplesize)*nstats));
  memset(REAL(prop_sample), 0, asInteger(samplesize)*nstats*sizeof(double));

  SEXP status;
  if(MHp) status = PROTECT(ScalarInteger(DISPATCH_SANSample(s,
                                                     REAL(invcov), REAL(tau), REAL(sample), REAL(prop_sample), asInteger(samplesize),
                                                     asInteger(nsteps),
                                                     nstats, INTEGER(statindices), noffsets, INTEGER(offsetindices), REAL(offsets),
                                                     asInteger(verbose))));
  else status = PROTECT(ScalarInteger(MCMC_MH_FAILED));

  const char *outn[] = {"status", "s", "s.prop", "state", ""};
  SEXP outl = PROTECT(mkNamed(VECSXP, outn));
  SET_VECTOR_ELT(outl, 0, status);
  SET_VECTOR_ELT(outl, 1, sample);
  SET_VECTOR_ELT(outl, 2, prop_sample);

  /* record new generated network to pass back to R */
  if(asInteger(status) == MCMC_OK){
    s->stats = REAL(sample) + (asInteger(samplesize)-1)*nstats;
    SET_VECTOR_ELT(outl, 3, DISPATCH_ErgmStateRSave(s));
  }

  DISPATCH_ErgmStateDestroy(s);
  PutRNGstate();  /* Disable RNG before returning */
  UNPROTECT(4);
  return outl;
}


/*********************
 void DISPATCH_SANSample

 Using the parameters contained in the array theta, obtain the
 network statistics for a sample of size samplesize.  nsteps is the
 initial number of Markov chain steps before sampling anything
 and interval is the number of MC steps between successive 
 networks in the sample.  Put all the sampled statistics into
 the networkstatistics array. 
*********************/
MCMCStatus DISPATCH_SANSample(DISPATCH_ErgmState *s,
                       double *invcov, double *tau, double *networkstatistics, double *prop_networkstatistics,
                       int samplesize, int nsteps,
                       int nstats,
                       int *statindices,
                       int noffsets,
                       int *offsetindices,
                       double *offsets,
                       int verbose){
  double *deltainvsig = R_calloc(nstats, double);

  int staken, tottaken, ptottaken;
  unsigned int interval = nsteps / samplesize; // Integer division: rounds down.
  unsigned int burnin = nsteps - (samplesize-1)*interval;

  if(DISPATCH_SANMetropolisHastings(s, invcov, tau, networkstatistics, prop_networkstatistics, burnin, &staken,
                                    nstats, statindices, noffsets, offsetindices, offsets, deltainvsig,
                             verbose)!=MCMC_OK)
    return MCMC_MH_FAILED;

  if (samplesize>1){
    staken = 0;
    tottaken = 0;
    ptottaken = 0;
    
    /* Now sample networks */
    for (unsigned int i=1; i < samplesize; i++){
      /* Set current vector of stats equal to previous vector */
      Rboolean finished = TRUE;
      for (unsigned int j=0; j<nstats; j++){
        if((networkstatistics[j+nstats] = networkstatistics[j])!=0) finished = FALSE;
      }
      if(finished){
	if(verbose) Rprintf("Exact match found.\n");
	break;
      }

      networkstatistics += nstats;
      prop_networkstatistics += nstats;
      /* This then adds the change statistics to these values */
      
      if(DISPATCH_SANMetropolisHastings(s, invcov, tau, networkstatistics, prop_networkstatistics,
                                 interval, &staken, nstats, statindices, noffsets, offsetindices, offsets,
                                        deltainvsig, verbose)!=MCMC_OK)
	return MCMC_MH_FAILED;
      tottaken += staken;
      if (verbose){
        if( ((3*i) % samplesize)==0 && samplesize > 500){
        Rprintf("Sampled %d from SAN Metropolis-Hastings\n", i);}
      }
      
      if( ((3*i) % samplesize)==0 && tottaken == ptottaken){
        ptottaken = tottaken; 
        Rprintf("Warning:  SAN Metropolis-Hastings algorithm has accepted only "
        "%d steps out of a possible %d\n",  ptottaken-tottaken, i); 
      }

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
	  Rprintf("SAN Metropolis-Hastings accepted %7.3f%% of %lld proposed steps.\n",
	    tottaken*100.0/(1.0*interval*samplesize), (long long) interval*samplesize); 
    }
  }else{
    if (verbose){
      Rprintf("SAN Metropolis-Hastings accepted %7.3f%% of %d proposed steps.\n",
	      staken*100.0/(1.0*nsteps), nsteps); 
    }
  }
  return MCMC_OK;
}

/*********************
MCMCStatus DISPATCH_SANMetropolisHastings

 In this function, theta is a m->n_stats-vector just as in SANSample,
 but now networkstatistics is merely another m->n_stats-vector because
 this function merely iterates nsteps times through the Markov
 chain, keeping track of the cumulative change statistics along
 the way, then returns, leaving the updated change statistics in
 the networkstatistics vector.  In other words, this function 
 essentially generates a sample of size one
*********************/
MCMCStatus DISPATCH_SANMetropolisHastings(DISPATCH_ErgmState *s,
                                   double *invcov,
                                   double *tau, double *networkstatistics, double *prop_networkstatistics,
                                   int nsteps, int *staken,
                                   int nstats,
                                   int *statindices,
                                   int noffsets,
                                   int *offsetindices,
                                   double *offsets,
                                          double *deltainvsig,
                                   int verbose){
  DISPATCH_Network *nwp = s->nwp;
  DISPATCH_Model *m = s->m;
  DISPATCH_MHProposal *MHp = s->MHp;

  unsigned int taken=0, unsuccessful=0;
  
/*  if (verbose)
    Rprintf("Now proposing %d DISPATCH_MH steps... ", nsteps); */
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

    /* Always store the proposal for self-tuning. */
    for (unsigned int i = 0; i < nstats; i++){
      prop_networkstatistics[i] += m->workspace[statindices[i]];
    }

    if(verbose>=5){
      Rprintf("Changes: (");
      for(unsigned int i=0; i<nstats; i++)
	Rprintf(" %f ", m->workspace[statindices[i]]);
      Rprintf(")\n");
    }
    
    /* Calculate the change in the (s-t) %*% W %*% (s-t) due to the proposal. */
    double ip=0;
    for (unsigned int i=0; i<nstats; i++){
     deltainvsig[i]=0.0;
     for (unsigned int j=0; j<nstats; j++){
      deltainvsig[i]+=(m->workspace[statindices[j]])*invcov[i+(nstats)*j];
     }
     ip+=deltainvsig[i]*((m->workspace[statindices[i]])+2.0*networkstatistics[i]);
    }
    
    double offsetcontrib = 0;
    for(int i = 0; i < noffsets; i++){
        offsetcontrib += (m->workspace[offsetindices[i]])*offsets[i];
    }
    
    if(verbose>=5){
      Rprintf("log acceptance probability: %f\n", -ip/tau[0] + offsetcontrib);
    }
    
    /* if we accept the proposed network */
    if (tau[0]==0? ip - offsetcontrib <= 0 : ip/tau[0] - offsetcontrib <= -log(unif_rand()) ) { 
      if(verbose>=5){
	Rprintf("Accepted.\n");
      }

      /* Make proposed toggles (updating timestamps--i.e., for real this time) */
      for(unsigned int i=0; i < MHp->ntoggles; i++) PROP_COMMIT;
      /* record network statistics for posterity */
      Rboolean finished = TRUE;
      for (unsigned int i = 0; i < nstats; i++){
	if((networkstatistics[i] += m->workspace[statindices[i]])!=0) finished = FALSE;
      }
      
      taken++;

      if(finished) break;
    }else{
      if(verbose>=5){
	Rprintf("Rejected.\n");
      }
    }
  }

  *staken = taken;
  return MCMC_OK;
}
