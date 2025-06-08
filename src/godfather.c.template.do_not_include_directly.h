/*  File src/godfather.c.template.do_not_include_directly.h in package ergm,
 *  part of the Statnet suite of packages for network analysis,
 *  https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */

#include "ergm_util.h"

MCMCStatus ETYPE(Godfather)(ETYPE(ErgmState) *s, Edge n_changes, Vertex *tails, Vertex *heads, EWTTYPE *weights, double *stats){
  ETYPE(Network) *nwp = s->nwp;
  ETYPE(Model) *m = s->m;

  /* Doing this one change at a time saves a lot of changes... */
  for(Edge e=0; e<n_changes; e++){
    Vertex t = tails[e], h = heads[e];

    if(t==0){
      memcpy(stats+m->n_stats, stats, m->n_stats*sizeof(double));
      stats+=m->n_stats;
      continue;
    }

    EWTTYPE w, edgestate = GETWT(t,h);
    if(weights){
      w = weights[e];
      if(edgestate==w) continue;
    }

    ETYPE(EXEC_THROUGH_TERMS_INTO)(m, stats, {
	if(mtp->c_func){
	  ZERO_ALL_CHANGESTATS();
	  (*(mtp->c_func))(t, h, IFEWT(w,)
			   mtp, nwp, edgestate);  /* Call c_??? function */
	}else if(mtp->d_func){
	  (*(mtp->d_func))(1, &t, &h, IFEWT(&w,)
			   mtp, nwp);  /* Call d_??? function */
	}
        addonto(dstats, mtp->dstats, N_CHANGE_STATS);
      });

    /* Update network */
    IFELSEEWT(ETYPE(SetEdge)(t, h, w, nwp),               \
                 ToggleKnownEdge(t, h, nwp, edgestate));
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
SEXP ETYPE(Godfather_wrapper)(SEXP stateR,
                         // Godfather settings
                         SEXP changetails, SEXP changeheads, SEXP changeweights,
                         SEXP end_network,
                         SEXP verbose){
  GetRNGstate();  /* R function enabling uniform RNG */
  ETYPE(ErgmState) *s = ETYPE(ErgmStateInit)(stateR, ERGM_STATE_NO_INIT_PROP);
  ETYPE(Model) *m = s->m;

  /* (# 0-sentinels) + 1 is the number of output rows. */
  unsigned int nstatrows = 1;
  for(int *ct = INTEGER(changetails), *cte = ct+length(changetails); ct < cte ; ct++) if(*ct==0) nstatrows++;

  SEXP stats = PROTECT(allocVector(REALSXP, m->n_stats*nstatrows));
  memcpy(REAL(stats), s->stats, m->n_stats*sizeof(double));

  SEXP status = PROTECT(ScalarInteger(ETYPE(Godfather)(s, length(changetails), (Vertex*)INTEGER(changetails), (Vertex*)INTEGER(changeheads), length(changeweights)==0? NULL : (EWTTYPE *) EWTRTYPE(changeweights), REAL(stats))));

  const char *outn[] = {"status", "s", "state", ""};
  SEXP outl = PROTECT(mkNamed(VECSXP, outn));
  SET_VECTOR_ELT(outl, 0, status);
  SET_VECTOR_ELT(outl, 1, stats);

  /* record new generated network to pass back to R */
  if(asInteger(status) == MCMC_OK && asInteger(end_network)){
    s->stats = REAL(stats) + (nstatrows-1)*m->n_stats;
    SET_VECTOR_ELT(outl, 2, ETYPE(ErgmStateRSave)(s));
  }

  ETYPE(ErgmStateDestroy)(s);
  PutRNGstate();  /* Disable RNG before returning */
  UNPROTECT(3);
  return outl;
}
