/*  File src/wtgodfather.c in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#include "wtMCMC.h"
#include "ergm_wtmodel.h"
#include "ergm_wtchangestat.h"
#include "ergm_wtstate.h"
#include "ergm_util.h"

MCMCStatus WtGodfather(WtErgmState *s, Edge n_changes, Vertex *tails, Vertex *heads, double *weights, double *stats){
  WtNetwork *nwp = s->nwp;
  WtModel *m = s->m;

  /* Doing this one change at a time saves a lot of changes... */
  for(Edge e=0; e<n_changes; e++){
    Vertex t = tails[e], h = heads[e];
    double w = weights[e], edgestate;

    if(t==0){
      memcpy(stats+m->n_stats, stats, m->n_stats*sizeof(double));
      stats+=m->n_stats;
      continue;
    }
    
    if((edgestate=GETWT(t,h))==w)
      continue;

    WtEXEC_THROUGH_TERMS_INTO(m, stats, {
	if(mtp->c_func){
	  ZERO_ALL_CHANGESTATS();
	  (*(mtp->c_func))(t, h, w,
			   mtp, nwp, edgestate);  /* Call c_??? function */
	}else if(mtp->d_func){
	  (*(mtp->d_func))(1, &t, &h, &w,
			   mtp, nwp);  /* Call d_??? function */
	}
        addonto(dstats, mtp->dstats, N_CHANGE_STATS);
      });


    /* Update network */
    WtSetEdge(t, h, w, nwp);
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
SEXP WtGodfather_wrapper(SEXP stateR,
                         // Godfather settings
                         SEXP changetails, SEXP changeheads, SEXP changeweights,
                         SEXP end_network,
                         SEXP verbose){
  GetRNGstate();  /* R function enabling uniform RNG */
  WtErgmState *s = WtErgmStateInit(stateR, ERGM_STATE_NO_INIT_PROP);
  WtModel *m = s->m;

  /* (# 0-sentinels) + 1 is the number of output rows. */
  unsigned int nstatrows = 1;
  for(int *ct = INTEGER(changetails), *cte = ct+length(changetails); ct < cte ; ct++) if(*ct==0) nstatrows++;

  SEXP stats = PROTECT(allocVector(REALSXP, m->n_stats*nstatrows));
  memcpy(REAL(stats), s->stats, m->n_stats*sizeof(double));

  SEXP status = PROTECT(ScalarInteger(WtGodfather(s, length(changetails), (Vertex*)INTEGER(changetails), (Vertex*)INTEGER(changeheads), REAL(changeweights), REAL(stats))));

  const char *outn[] = {"status", "s", "state", ""};
  SEXP outl = PROTECT(mkNamed(VECSXP, outn));
  SET_VECTOR_ELT(outl, 0, status);
  SET_VECTOR_ELT(outl, 1, stats);

  /* record new generated network to pass back to R */
  if(asInteger(status) == MCMC_OK && asInteger(end_network)){
    s->stats = REAL(stats) + (nstatrows-1)*m->n_stats;
    SET_VECTOR_ELT(outl, 2, WtErgmStateRSave(s));
  }
  
  WtErgmStateDestroy(s);
  PutRNGstate();  /* Disable RNG before returning */
  UNPROTECT(3);
  return outl;
}
