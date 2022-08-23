/*  File src/godfather.c in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2022 Statnet Commons
 */
#include "MCMC.h"
#include "ergm_model.h"
#include "ergm_changestat.h"
#include "ergm_state.h"
#include "ergm_util.h"

MCMCStatus Godfather(ErgmState *s, Edge n_changes, Vertex *tails, Vertex *heads, int *weights, double *stats){
  Network *nwp = s->nwp;
  Model *m = s->m;

  /* Doing this one change at a time saves a lot of changes... */
  for(Edge e=0; e<n_changes; e++){
    Vertex t=TAIL(e), h=HEAD(e);

    if(t==0){
      memcpy(stats+m->n_stats, stats, m->n_stats*sizeof(double));
      stats+=m->n_stats;
      continue;
    }
    
    Rboolean edgestate = IS_OUTEDGE(t,h);
    if(weights){
      if(edgestate==weights[e])
	continue;
    }

    EXEC_THROUGH_TERMS_INTO(m, stats, {
	if(mtp->c_func){
	  ZERO_ALL_CHANGESTATS();
	  (*(mtp->c_func))(t, h,
			   mtp, nwp, edgestate);  /* Call c_??? function */
	}else if(mtp->d_func){
	  (*(mtp->d_func))(1, &t, &h,
			   mtp, nwp);  /* Call d_??? function */
	}
        addonto(dstats, mtp->dstats, N_CHANGE_STATS);
      });

    /* Update network */
    ToggleKnownEdge(t, h, nwp, edgestate);
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
SEXP Godfather_wrapper(SEXP stateR,
                       // Godfather settings
                       SEXP changetails, SEXP changeheads, SEXP changeweights,
                       SEXP end_network,
                       SEXP verbose){
  GetRNGstate();  /* R function enabling uniform RNG */
  ErgmState *s = ErgmStateInit(stateR, ERGM_STATE_NO_INIT_PROP);
  Model *m = s->m;

  /* (# 0-sentinels) + 1 is the number of output rows. */
  unsigned int nstatrows = 1;
  for(int *ct = INTEGER(changetails), *cte = ct+length(changetails); ct < cte ; ct++) if(*ct==0) nstatrows++;

  SEXP stats = PROTECT(allocVector(REALSXP, m->n_stats*nstatrows));
  memcpy(REAL(stats), s->stats, m->n_stats*sizeof(double));

  SEXP status = PROTECT(ScalarInteger(Godfather(s, length(changetails), (Vertex*)INTEGER(changetails), (Vertex*)INTEGER(changeheads),
                                                length(changeweights)==0? NULL : INTEGER(changeweights), REAL(stats))));

  const char *outn[] = {"status", "s", "state", ""};
  SEXP outl = PROTECT(mkNamed(VECSXP, outn));
  SET_VECTOR_ELT(outl, 0, status);
  SET_VECTOR_ELT(outl, 1, stats);

  /* record new generated network to pass back to R */
  if(asInteger(status) == MCMC_OK && asInteger(end_network)){
    s->stats = REAL(stats) + (nstatrows-1)*m->n_stats;
    SET_VECTOR_ELT(outl, 2, ErgmStateRSave(s));
  }

  ErgmStateDestroy(s);
  PutRNGstate();  /* Disable RNG before returning */
  UNPROTECT(3);
  return outl;
}
