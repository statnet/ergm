/*  File src/CD.c.template.do_not_include_directly.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2024 Statnet Commons
 */
#include "ergm_util.h"
/*****************
 Note on undirected networks:  For j<k, edge {j,k} should be stored
 as (j,k) rather than (k,j).  In other words, only directed networks
 should have (k,j) with k>j.
*****************/

/*****************
 void DISPATCH_CD_wrapper

 Wrapper for a call from R.

 and don't forget that tail -> head
*****************/
SEXP DISPATCH_CD_wrapper(SEXP stateR,
                  // MCMC settings
                  SEXP eta, SEXP samplesize, 
                  SEXP CDparams,
                  SEXP verbose){
  GetRNGstate();  /* R function enabling uniform RNG */
  DISPATCH_ErgmState *s = DISPATCH_ErgmStateInit(stateR, 0);

  DISPATCH_Model *m = s->m;
  DISPATCH_MHProposal *MHp = s->MHp;

  CD_UNDOS_ALLOC;
  double *extraworkspace = R_calloc(m->n_stats, double);

  SEXP sample = PROTECT(allocVector(REALSXP, asInteger(samplesize)*m->n_stats));
  memset(REAL(sample), 0, asInteger(samplesize)*m->n_stats*sizeof(double));

  SEXP status;
  if(MHp) status = PROTECT(ScalarInteger(DISPATCH_CDSample(s,
                                                    REAL(eta), REAL(sample), asInteger(samplesize), INTEGER(CDparams), CD_UNDOS_PASS, extraworkspace,
                                                    asInteger(verbose))));
  else status = PROTECT(ScalarInteger(MCMC_MH_FAILED));

  const char *outn[] = {"status", "s", ""};
  SEXP outl = PROTECT(mkNamed(VECSXP, outn));
  SET_VECTOR_ELT(outl, 0, status);
  SET_VECTOR_ELT(outl, 1, sample);

  DISPATCH_ErgmStateDestroy(s);  
  PutRNGstate();  /* Disable RNG before returning */
  UNPROTECT(3);
  return outl;
}


/*********************
 void DISPATCH_CDSample

 Using the parameters contained in the array eta, obtain the
 network statistics for a sample of size samplesize.  burnin is the
 initial number of Markov chain steps before sampling anything
 and interval is the number of MC steps between successive 
 networks in the sample.  Put all the sampled statistics into
 the networkstatistics array. 
*********************/
MCMCStatus DISPATCH_CDSample(DISPATCH_ErgmState *s,
                        double *eta, double *networkstatistics, 
			int samplesize, int *CDparams,
                      CD_UNDOS_RECEIVE, double *extraworkspace,
                        int verbose){
  DISPATCH_Model *m = s->m;

  int staken=0;
    
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

  /* Now sample networks */
  unsigned int i=0, sattempted=0;
  while(i<samplesize){
    
    if(DISPATCH_CDStep(s, eta, networkstatistics, CDparams, &staken, CD_UNDOS_PASS, extraworkspace,
		verbose)!=MCMC_OK)
      return MCMC_MH_FAILED;

    R_CheckUserInterruptEvery(16L, i);
#ifdef Win32
    if( ((100*i) % samplesize)==0 && samplesize > 500){
      R_FlushConsole();
      R_ProcessEvents();
    }
#endif

      networkstatistics += m->n_stats;
      i++;

    sattempted++;
  }

  if (verbose){
    Rprintf("Sampler accepted %7.3f%% of %lld proposed steps.\n",
	    staken*100.0/(1.0*sattempted*CDparams[0]), (long long) sattempted*CDparams[0]); 
  }
  
  return MCMC_OK;
}

/*********************
 void MetropolisHastings

 In this function, eta is a m->n_stats-vector just as in DISPATCH_CDSample,
 but now networkstatistics is merely another m->n_stats-vector because
 this function merely iterates nsteps=CDparams[0] times through the Markov
 chain, keeping track of the cumulative change statistics along
 the way, then returns, leaving the updated change statistics in
 the networkstatistics vector.  In other words, this function 
 essentially generates a sample of size one
*********************/
MCMCStatus DISPATCH_CDStep(DISPATCH_ErgmState *s,
                      double *eta, double *networkstatistics,
                      int *CDparams, int *staken,
                      CD_UNDOS_RECEIVE, double *extraworkspace,
                      int verbose){

  DISPATCH_Network *nwp = s->nwp;
  DISPATCH_Model *m = s->m;
  DISPATCH_MHProposal *MHp = s->MHp;
  
  unsigned int unsuccessful=0, ntoggled=0;

  for(unsigned int step=0; step<CDparams[0]; step++){
    unsigned int mtoggled=0;
    memset(extraworkspace, 0, m->n_stats*sizeof(double));
    double cumlr = 0;
    
    for(unsigned int mult=0; mult<CDparams[1]; mult++){
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
	  if(unsuccessful>*staken*MH_QUIT_UNSUCCESSFUL){
	    Rprintf("Too many MH MHProposal function failures.\n");
	    return MCMC_MH_FAILED;
	  }
	  continue;
	  
	case MH_CONSTRAINT:
	  cumlr = MHp->logratio = -INFINITY; // Force rejection of proposal.
	  goto REJECT;
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

      // Add them to the cumulative changes.
      addonto(extraworkspace, m->workspace, m->n_stats);
      
      if(verbose>=5){
	Rprintf("Changes: (");
	for(unsigned int i=0; i<m->n_stats; i++){
	  Rprintf(" %f ", m->workspace[i]);
	}
	Rprintf(")\n");
      }

      if(mult<CDparams[1]-1){
	/* Make proposed toggles provisionally. */
	for(unsigned int i=0; i < MHp->ntoggles; i++){
	  CD_PROP_TOGGLE_PROVISIONAL;
          mtoggled++;
	}
      }

      // Accumulate the log acceptance ratio.
      cumlr += MHp->logratio;
    } // mult

    
    if(verbose>=5){
      Rprintf("Cumulative changes: (");
      for(unsigned int i=0; i<m->n_stats; i++)
	Rprintf(" %f ", extraworkspace[i]);
      Rprintf(")\n");
    }
    
    /* Calculate inner (dot) product */
    double ip = dotprod(eta, extraworkspace, m->n_stats);

    /* The logic is to set cutoff = ip+logratio ,
       then let the MH probability equal min{exp(cutoff), 1.0}.
       But we'll do it in log space instead.  */
    double cutoff = ip + cumlr;

    if(verbose>=5){
      Rprintf("log acceptance probability: %f + %f = %f\n", ip, cumlr, cutoff);
    }
    
    /* if we accept the proposed network */
    if (cutoff >= 0.0 || logf(unif_rand()) < cutoff) { 
      if(verbose>=5){
	Rprintf("Accepted.\n");
      }
      (*staken)++; 

      if(step<CDparams[0]-1){
	/* Make the remaining proposed toggles (which we did not make provisionally) */
	/* Then, make the changes. */
	for(unsigned int i=0; i < MHp->ntoggles; i++){
	  CD_PROP_TOGGLE_PROVISIONAL;
	}
      }

      /* record network statistics for posterity */
      addonto(networkstatistics, extraworkspace, m->n_stats);

    }else{
    REJECT:
      if(verbose>=5){
	Rprintf("Rejected.\n");
      }
      // Undo the provisional toggles (the last mtoggled ones)
      for(unsigned int i=0; i < mtoggled; i++){
	ntoggled--;
	CD_PROP_UNDO_TOGGLE(ntoggled);
      }
    }
  } // step
  
  /* Undo toggles. */
  while(ntoggled){
    ntoggled--;
    CD_PROP_UNDO_TOGGLE(ntoggled)
  }

  return MCMC_OK;
}

