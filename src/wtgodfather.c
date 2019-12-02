/*  File src/wtgodfather.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2019 Statnet Commons
 */
#include "wtMCMC.h"
#include "ergm_wtmodel.h"
#include "ergm_wtchangestat.h"
#include "ergm_wtstate.h"
#include "ergm_util.h"

MCMCStatus WtGodfather(ErgmWtState *s, Edge n_changes, Vertex *tails, Vertex *heads, double *weights, double *stats){
  WtNetwork *nwp = s->nwp;
  WtModel *m = s->m;

  stats+=m->n_stats;

  /* Doing this one change at a time saves a lot of changes... */
  for(Edge e=0; e<n_changes; e++){
    Vertex t = tails[e], h = heads[e];
    double w = weights[e], edgeweight;

    if(t==0){
      stats+=m->n_stats;
      continue;
    }
    
    if((edgeweight=GETWT(t,h))==w)
      continue;

    WtEXEC_THROUGH_TERMS_INTO(m, stats, {
	if(mtp->c_func){
	  ZERO_ALL_CHANGESTATS();
	  (*(mtp->c_func))(t, h, w,
			   mtp, nwp, edgeweight);  /* Call c_??? function */
	}else if(mtp->d_func){
	  (*(mtp->d_func))(1, &t, &h, &w,
			   mtp, nwp);  /* Call d_??? function */
	}
	for(unsigned int k=0; k<N_CHANGE_STATS; k++){
	    dstats[k] += mtp->dstats[k];
	}
      });


    /* Update storage and network */    
    WtUPDATE_STORAGE_SET(t, h, w, nwp, m, NULL, edgeweight);
  }

  return MCMC_OK;
}

/*****************
 void Godfather_wrapper

 ...we'll make them an offer (of changes) they can't refuse.
 This function takes a list of changes, each with a time stamp,
 then produces a matrix of changestats (with one row for each unique
 time stamp value) that result from performing all the changes at
 each time step.  For instance, one might use this function to 
 find the changestats that result from starting from an empty network
 and then adding all of the edges to make up an observed network of interest.
*****************/
SEXP WtGodfather_wrapper(// Network settings
                         SEXP dn, SEXP dflag, SEXP bipartite,
                         // Model settings
                         SEXP nterms, SEXP funnames,
                         SEXP sonames,
                         // Numeric inputs
                         SEXP inputs,
                         // Network state
                         SEXP nedges,
                         SEXP tails, SEXP heads, SEXP weights,
                         // Godfather settings
                         SEXP nsteps,
                         SEXP changetails, SEXP changeheads, SEXP changeweights,
                         SEXP end_network,
                         SEXP verbose){
  GetRNGstate();  /* R function enabling uniform RNG */
  ErgmWtState *s = ErgmWtStateInit(// Network settings
                               asInteger(dn), asInteger(dflag), asInteger(bipartite),
                               // Model settings
                               asInteger(nterms), FIRSTCHAR(funnames), FIRSTCHAR(sonames), FALSE,
                               // Proposal settings
                               NO_WTMHPROPOSAL,
                               // Numeric inputs
                               REAL(inputs),
                               // Network state
                               asInteger(nedges), (Vertex*) INTEGER(tails), (Vertex*) INTEGER(heads), REAL(weights),
                               NO_LASTTOGGLE);
  WtNetwork *nwp = s->nwp;
  WtModel *m = s->m;

  SEXP stats = PROTECT(allocVector(REALSXP, m->n_stats*(1+asInteger(nsteps))));
  memset(REAL(stats), 0, m->n_stats*(1+asInteger(nsteps))*sizeof(double));

  SEXP status = PROTECT(ScalarInteger(WtGodfather(s, length(changetails), (Vertex*)INTEGER(changetails), (Vertex*)INTEGER(changeheads), REAL(changeweights), REAL(stats))));

  const char *outn[] = {"status", "s", WTNWSTATE_NAMES, ""};
  SEXP outl = PROTECT(mkNamed(VECSXP, outn));
  SET_VECTOR_ELT(outl, 0, status);
  SET_VECTOR_ELT(outl, 1, stats);

  /* record new generated network to pass back to R */
  if(asInteger(status) == MCMC_OK && asInteger(end_network)){
    WTNWSTATE_SAVE_INTO_RLIST(nwp, outl, 2);
  }
  
  ErgmWtStateDestroy(s);
  PutRNGstate();  /* Disable RNG before returning */
  UNPROTECT(3);
  return outl;
}
