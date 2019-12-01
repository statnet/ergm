/*  File src/wtCD.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2019 Statnet Commons
 */
#include "wtCD.h"
#include "ergm_util.h"

/*****************
 Note on undirected networks:  For j<k, edge {j,k} should be stored
 as (j,k) rather than (k,j).  In other words, only directed networks
 should have (k,j) with k>j.
*****************/

/*****************
 void WtCD_wrapper

 Wrapper for a call from R.

 and don't forget that tail -> head
*****************/
SEXP WtCD_wrapper(// Network settings
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
                  SEXP tails, SEXP heads, SEXP weights,
                  // MCMC settings
                  SEXP eta, SEXP samplesize, 
                  SEXP CDparams,
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

  WtModel *m = s->m;
  WtMHProposal *MHp = s->MHp;

  Vertex *undotail = Calloc(MHp->ntoggles * INTEGER(CDparams)[0] * INTEGER(CDparams)[1], Vertex);
  Vertex *undohead = Calloc(MHp->ntoggles * INTEGER(CDparams)[0] * INTEGER(CDparams)[1], Vertex);
  double *undoweight = Calloc(MHp->ntoggles * INTEGER(CDparams)[0] * INTEGER(CDparams)[1], double);
  double *extraworkspace = Calloc(m->n_stats, double);

  SEXP sample = PROTECT(allocVector(REALSXP, asInteger(samplesize)*m->n_stats));
  memset(REAL(sample), 0, asInteger(samplesize)*m->n_stats*sizeof(double));

  SEXP status;
  if(MHp) status = PROTECT(ScalarInteger(WtCDSample(s,
                                                    REAL(eta), REAL(sample), asInteger(samplesize), INTEGER(CDparams), undotail, undohead, undoweight, extraworkspace,
                                                    asInteger(verbose))));
  else status = PROTECT(ScalarInteger(WtMCMC_MH_FAILED));

  SEXP outl = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(outl, 0, status);
  SET_VECTOR_ELT(outl, 1, sample);

  Free(undotail);
  Free(undohead);
  Free(undoweight);
  Free(extraworkspace);

  ErgmWtStateDestroy(s);  
  PutRNGstate();  /* Disable RNG before returning */
  UNPROTECT(3);
  return outl;
}


/*********************
 void WtCDSample

 Using the parameters contained in the array eta, obtain the
 network statistics for a sample of size samplesize.  burnin is the
 initial number of Markov chain steps before sampling anything
 and interval is the number of MC steps between successive 
 networks in the sample.  Put all the sampled statistics into
 the networkstatistics array. 
*********************/
WtMCMCStatus WtCDSample(ErgmWtState *s,
                        double *eta, double *networkstatistics, 
			int samplesize, int *CDparams,
                        Vertex *undotail, Vertex *undohead, double *undoweight, double *extraworkspace,
                        int verbose){
  WtModel *m = s->m;

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
    
    if(WtCDStep(s, eta, networkstatistics, CDparams, &staken, undotail, undohead, undoweight, extraworkspace,
		verbose)!=WtMCMC_OK)
      return WtMCMC_MH_FAILED;
    
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

  return WtMCMC_OK;
}

/*********************
 void MetropolisHastings

 In this function, eta is a m->n_stats-vector just as in WtCDSample,
 but now networkstatistics is merely another m->n_stats-vector because
 this function merely iterates nsteps=CDparams[0] times through the Markov
 chain, keeping track of the cumulative change statistics along
 the way, then returns, leaving the updated change statistics in
 the networkstatistics vector.  In other words, this function 
 essentially generates a sample of size one
*********************/
WtMCMCStatus WtCDStep(ErgmWtState *s,
                      double *eta, double *networkstatistics,
                      int *CDparams, int *staken,
                      Vertex *undotail, Vertex *undohead, double *undoweight, double *extraworkspace,
                      int verbose){

  WtNetwork *nwp = s->nwp;
  WtModel *m = s->m;
  WtMHProposal *MHp = s->MHp;

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
	  return WtMCMC_MH_FAILED;
	  
	case MH_UNSUCCESSFUL:
	  warning("MH MHProposal function failed to find a valid proposal.");
	  unsuccessful++;
	  if(unsuccessful>*staken*MH_QUIT_UNSUCCESSFUL){
	    Rprintf("Too many MH MHProposal function failures.\n");
	    return WtMCMC_MH_FAILED;
	  }
	  continue;
	  
	case MH_CONSTRAINT:
	  cumlr = MHp->logratio = -INFINITY; // Force rejection of proposal.
	  goto REJECT;
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

      // Add them to the cumulative changes.
      for(unsigned int i=0; i<m->n_stats; i++)
	extraworkspace[i] += m->workspace[i];
      
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
	  Vertex t=MHp->toggletail[i], h=MHp->togglehead[i];
	  double w=MHp->toggleweight[i];
	  undotail[ntoggled]=t;
	  undohead[ntoggled]=h;
	  undoweight[ntoggled]=WtGetEdge(MHp->toggletail[i], MHp->togglehead[i], nwp);
	  ntoggled++;
	  mtoggled++;

	  WtUPDATE_STORAGE_SET(t, h, w, nwp, m, MHp, undoweight[ntoggled]);
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
	  Vertex t=MHp->toggletail[i], h=MHp->togglehead[i];
	  double w=MHp->toggleweight[i];
	  undotail[ntoggled]=t;
	  undohead[ntoggled]=h;
	  undoweight[ntoggled]=WtGetEdge(MHp->toggletail[i], MHp->togglehead[i], nwp);
	  ntoggled++;

	  WtUPDATE_STORAGE_SET(t, h, w, nwp, m, MHp, undoweight[ntoggled]);
	}
      }

      /* record network statistics for posterity */
      for (unsigned int i = 0; i < m->n_stats; i++){
	networkstatistics[i] += extraworkspace[i];
      }

    }else{
    REJECT:
      if(verbose>=5){
	Rprintf("Rejected.\n");
      }
      // Undo the provisional toggles (the last mtoggled ones)
      for(unsigned int i=0; i < mtoggled; i++){
	ntoggled--;
	Vertex t = undotail[ntoggled], h = undohead[ntoggled];
	double w = undoweight[ntoggled];

	/* FIXME: This should be done in one call, but it's very easy
	   to make a fencepost error here. */
	WtGET_EDGE_UPDATE_STORAGE_SET(t, h, w, nwp, m, MHp);
      }
    }
  } // step
  
  /* Undo toggles. */
  for(unsigned int i=0; i < ntoggled; i++){
    Vertex t = undotail[i], h = undohead[i];
    double w = undoweight[i];

    /* FIXME: This should be done in one call, but it's very easy
       to make a fencepost error here. */
    WtGET_EDGE_UPDATE_STORAGE_SET(t, h, w, nwp, m, MHp);
  }
  
  return WtMCMC_OK;
}

