/*  File src/wtnetstats.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2019 Statnet Commons
 */
#include "wtnetstats.h"
#include "ergm_omp.h"
#include "ergm_util.h"
/*****************
 void network_stats_wrapper

 Wrapper for a call from R.  Return the statistics when of the observed graph.
*****************/

SEXP wt_network_stats_wrapper(ARGS_WTSTATE){
  GetRNGstate();  /* R function enabling uniform RNG */
  WtErgmState *s = WtErgmStateInit(YES_WTSTATE_EMPTY_NO_INIT_S);

  WtModel *m = s->m;

  SEXP stats = PROTECT(allocVector(REALSXP, m->n_stats));

  SEXP elR = getListElement(stateR, "el");
  WtSummStats(s,
              length(VECTOR_ELT(elR, 0)),
              (Vertex*) INTEGER(VECTOR_ELT(elR, 0)),
              (Vertex*) INTEGER(VECTOR_ELT(elR, 1)),
              REAL(VECTOR_ELT(elR, 2)),
              REAL(stats));
  
  PutRNGstate();
  UNPROTECT(1);
  return stats;
}

void WtEmptyNetworkStats(WtErgmState *s, Rboolean skip_s, double *stats){
  WtModel *m = s->m;

  WtEXEC_THROUGH_TERMS_INTO(m, stats, {
      if(!skip_s || mtp->s_func==NULL){
        if(mtp->emptynwstats)
          memcpy(dstats, mtp->emptynwstats, mtp->nstats*sizeof(double));
      }});
}

/****************
 void WtSummStats Computes summary statistics for a network. Must be
 passed an empty network.
*****************/
void WtSummStats(WtErgmState *s, Edge n_edges, Vertex *tails, Vertex *heads, double *weights, double *stats){
  WtNetwork *nwp = s->nwp;
  WtModel *m = s->m;

  memset(stats, 0, m->n_stats*sizeof(double));

  WtEmptyNetworkStats(s, TRUE, stats);
  WtZStats(nwp, m);
  addonto(stats, m->workspace, m->n_stats);
  

  WtDetShuffleEdges(tails,heads,weights,n_edges); /* Shuffle edgelist. */
  
  Edge ntoggles = n_edges; // So that we can use the macros

  /* Calculate statistics for terms that don't have c_functions or s_functions.  */
  WtEXEC_THROUGH_TERMS_INTO(m, stats, {
      if(mtp->s_func==NULL && mtp->c_func==NULL && mtp->d_func){
	(*(mtp->d_func))(ntoggles, tails, heads, weights,
			 mtp, nwp);  /* Call d_??? function */
	addonto(dstats, mtp->dstats, N_CHANGE_STATS);
      }
    });

  /* Calculate statistics for terms that have c_functions but not s_functions.  */
  FOR_EACH_TOGGLE{
    GETNEWTOGGLEINFO();
    
    ergm_PARALLEL_FOR_LIMIT(m->n_terms)
    WtEXEC_THROUGH_TERMS_INTO(m, stats, {
	if(mtp->s_func==NULL && mtp->c_func){
	  ZERO_ALL_CHANGESTATS();
	  (*(mtp->c_func))(TAIL, HEAD, NEWWT,
			   mtp, nwp, 0);  /* Call c_??? function */
	    addonto(dstats, mtp->dstats, N_CHANGE_STATS);
	}
      });
    
    /* Update storage and network */    
    WtUPDATE_STORAGE_COND(TAIL, HEAD, NEWWT, nwp, m, NULL, 0, mtp->s_func==NULL && mtp->d_func==NULL);
    SETWT(TAIL, HEAD, NEWWT);
  }
  
  /* Calculate statistics for terms have s_functions  */
  WtEXEC_THROUGH_TERMS_INTO(m, stats, {
      if(mtp->s_func){
	ZERO_ALL_CHANGESTATS();
	(*(mtp->s_func))(mtp, nwp);  /* Call d_??? function */
	for(unsigned int k=0; k<N_CHANGE_STATS; k++){
	  dstats[k] = mtp->dstats[k]; // Overwrite, not accumulate.
	}
      }
    });
}
