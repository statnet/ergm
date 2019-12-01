/*  File src/wtSAN.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2019 Statnet Commons
 */
#include "wtSAN.h"
#include "ergm_util.h"

/*****************
 Note on undirected networks:  For j<k, edge {j,k} should be stored
 as (j,k) rather than (k,j).  In other words, only directed networks
 should have (k,j) with k>j.
*****************/

/*****************
 void WtSAN_wrapper

 Wrapper for a call from R.
*****************/

SEXP WtSAN_wrapper(// Network settings
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
                   SEXP tau, SEXP stats,
                   SEXP samplesize, SEXP nsteps,
                   SEXP invcov,
                   SEXP maxedges,
                   SEXP nstats,
                   SEXP statindices,
                   SEXP noffsets,
                   SEXP offsetindices,
                   SEXP offsets,
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
  WtMHProposal *MHp = s->MHp;

  SEXP sample = PROTECT(allocVector(REALSXP, asInteger(samplesize)*asInteger(nstats)));
  memcpy(REAL(sample), REAL(stats), asInteger(nstats)*sizeof(double));
  SEXP prop_sample = PROTECT(allocVector(REALSXP, asInteger(samplesize)*asInteger(nstats)));
  memset(REAL(prop_sample), 0, asInteger(samplesize)*asInteger(nstats)*sizeof(double));

  SEXP status;
  if(MHp) status = PROTECT(ScalarInteger(WtSANSample(s,
                                                     REAL(invcov), REAL(tau), REAL(sample), REAL(prop_sample), asInteger(samplesize),
                                                     asInteger(nsteps),
                                                     asInteger(maxedges),
                                                     asInteger(nstats), INTEGER(statindices), asInteger(noffsets), INTEGER(offsetindices), REAL(offsets),
                                                     asInteger(verbose))));
  else status = PROTECT(ScalarInteger(WtMCMC_MH_FAILED));

  const char *outn[] = {"status", "s", "s.prop", "newnwtails", "newnwheads", "newnwweights", ""};
  SEXP outl = PROTECT(mkNamed(VECSXP, outn));
  SET_VECTOR_ELT(outl, 0, status);
  SET_VECTOR_ELT(outl, 1, sample);
  SET_VECTOR_ELT(outl, 2, prop_sample);

  /* record new generated network to pass back to R */
  if(asInteger(status) == WtMCMC_OK && asInteger(maxedges)>0){
    SEXP newnetworktails = PROTECT(allocVector(INTSXP, EDGECOUNT(nwp)+1));
    SEXP newnetworkheads = PROTECT(allocVector(INTSXP, EDGECOUNT(nwp)+1));
    SEXP newnetworkweights = PROTECT(allocVector(REALSXP, EDGECOUNT(nwp)+1));

    INTEGER(newnetworktails)[0]=INTEGER(newnetworkheads)[0]=REAL(newnetworkweights)[0]=
      WtEdgeTree2EdgeList((Vertex*)INTEGER(newnetworktails)+1,
                          (Vertex*)INTEGER(newnetworkheads)+1,
                          REAL(newnetworkweights)+1,
                          nwp,asInteger(maxedges)-1);

    SET_VECTOR_ELT(outl, 3, newnetworktails);
    SET_VECTOR_ELT(outl, 4, newnetworkheads);
    SET_VECTOR_ELT(outl, 5, newnetworkweights);
    UNPROTECT(3);
  }

  ErgmWtStateDestroy(s);
  PutRNGstate();  /* Disable RNG before returning */
  UNPROTECT(4);
  return outl;
}


/*********************
 void WtSANSample

 Using the parameters contained in the array theta, obtain the
 network statistics for a sample of size samplesize.  nsteps is the
 initial number of Markov chain steps before sampling anything
 and interval is the number of MC steps between successive 
 networks in the sample.  Put all the sampled statistics into
 the networkstatistics array. 
*********************/
WtMCMCStatus WtSANSample(ErgmWtState *s,
                         double *invcov, double *tau, double *networkstatistics, double *prop_networkstatistics,
                         int samplesize, int nsteps, 
                         int nmax,
                         int nstats,
                         int *statindices,
                         int noffsets,
                         int *offsetindices,
                         double *offsets,
                         int verbose){
  WtNetwork *nwp = s->nwp;

  int staken, tottaken, ptottaken;
    
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

  unsigned int interval = nsteps / samplesize; // Integer division: rounds down.
  unsigned int burnin = nsteps - (samplesize-1)*interval;
  
  /*********************
   Burn in step.  While we're at it, use nsteps statistics to 
   prepare covariance matrix for Mahalanobis distance calculations 
   in subsequent calls to M-H
   *********************/
  /*  Catch more edges than we can return */
  if(WtSANMetropolisHastings(s, invcov, tau, networkstatistics, prop_networkstatistics, burnin, &staken,
			     nstats, statindices, noffsets, offsetindices, offsets,
                             verbose)!=WtMCMC_OK)
    return WtMCMC_MH_FAILED;
  if(nmax!=0 && EDGECOUNT(nwp) >= nmax-1){
    return WtMCMC_TOO_MANY_EDGES;
  }

  if (samplesize>1){
    staken = 0;
    tottaken = 0;
    ptottaken = 0;
    
    /* Now sample networks */
    for (unsigned int i=1; i < samplesize; i++){
      /* Set current vector of stats equal to previous vector */
      Rboolean found = TRUE;
      for (unsigned int j=0; j<nstats; j++){
        if((networkstatistics[j+nstats] = networkstatistics[j])!=0) found = FALSE;
      }
      if(found){
	if(verbose) Rprintf("Exact match found.\n");
	break;
      }

      networkstatistics += nstats;
      prop_networkstatistics += nstats;
      /* This then adds the change statistics to these values */
      
      if(WtSANMetropolisHastings(s, invcov, tau, networkstatistics, prop_networkstatistics,
                                 interval, &staken, nstats, statindices, noffsets, offsetindices, offsets,
                                 verbose)!=WtMCMC_OK)
	return WtMCMC_MH_FAILED;
      if(nmax!=0 && EDGECOUNT(nwp) >= nmax-1){
	return WtMCMC_TOO_MANY_EDGES;
      }
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
  return WtMCMC_OK;
}

/*********************
MCMCStatus WtSANMetropolisHastings

 In this function, theta is a m->n_stats-vector just as in SANSample,
 but now networkstatistics is merely another m->n_stats-vector because
 this function merely iterates nsteps times through the Markov
 chain, keeping track of the cumulative change statistics along
 the way, then returns, leaving the updated change statistics in
 the networkstatistics vector.  In other words, this function 
 essentially generates a sample of size one
*********************/
WtMCMCStatus WtSANMetropolisHastings(ErgmWtState *s,
                                     double *invcov, 
                                     double *tau, double *networkstatistics, double *prop_networkstatistics,
                                     int nsteps, int *staken,
                                     int nstats,
                                     int *statindices,
                                     int noffsets,
                                     int *offsetindices,
                                     double *offsets,
                                     int verbose){
  WtNetwork *nwp = s->nwp;
  WtModel *m = s->m;
  WtMHProposal *MHp = s->MHp;

  unsigned int taken=0, unsuccessful=0;
  double *deltainvsig;
  deltainvsig = (double *)Calloc(nstats, double);
  
/*  if (verbose)
    Rprintf("Now proposing %d WtMH steps... ", nsteps); */
  for(unsigned int step=0; step < nsteps; step++) {
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
	if(unsuccessful>taken*MH_QUIT_UNSUCCESSFUL){
	  Rprintf("Too many MH MHProposal function failures.\n");
	  return WtMCMC_MH_FAILED;
	}
      case MH_CONSTRAINT:
	continue;
      }
    }
    
    if(verbose>=5){
      Rprintf("MHProposal: ");
      for(unsigned int i=0; i<MHp->ntoggles; i++)
	Rprintf(" (%d, %d)", MHp->toggletail[i], MHp->togglehead[i]);
      Rprintf("\n");
    }

    /* Calculate change statistics,
     remembering that tail -> head */
    WtChangeStats(MHp->ntoggles, MHp->toggletail, MHp->togglehead, MHp->toggleweight, nwp, m);

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
      Rprintf("log acceptance probability: %f\n", ip - offsetcontrib);
    }
    
    /* if we accept the proposed network */
    if (tau[0]==0? ip - offsetcontrib <= 0 : ip/tau[0] - offsetcontrib <= -log(unif_rand()) ) { 
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
      Rboolean found = TRUE;
      for (unsigned int i = 0; i < nstats; i++){
	if((networkstatistics[i] += m->workspace[statindices[i]])!=0) found=FALSE;
      }
      
      taken++;

      if(found)	break;
    }else{
      if(verbose>=5){
	Rprintf("Rejected.\n");
      }
    }
  }

  Free(deltainvsig);

  *staken = taken;
  return WtMCMC_OK;
}
