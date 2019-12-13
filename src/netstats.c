/*  File src/netstats.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2019 Statnet Commons
 */
#include "netstats.h"
#include "ergm_omp.h"
#include "ergm_util.h"
/*****************
 void network_stats_wrapper

 Wrapper for a call from R.  Return the change in the statistics when
 we go from an empty graph to the observed graph.  If the empty graph
 has true global values equal to zero for all statistics, then this
 change gives the true global values for the observed graph.
*****************/

SEXP network_stats_wrapper(ARGS_STATE,
                           ARGS_LASTTOGGLE){
  GetRNGstate();  /* R function enabling uniform RNG */
  ErgmState *s = ErgmStateInit(YES_STATE_EMPTY_NO_INIT_S,
                               YES_LASTTOGGLE);

  Model *m = s->m;

  SEXP stats = PROTECT(allocVector(REALSXP, m->n_stats));

  if(s->stats) memcpy(REAL(stats), s->stats, m->n_stats*sizeof(double));
  else memset(REAL(stats), 0, m->n_stats*sizeof(double));

  /* Compute the change statistics and copy them to stats for return
     to R.  Note that stats already has the statistics of an empty
     network, so d_??? statistics will add on to them, while s_???
     statistics will simply overwrite them. */
  SEXP elR = getListElement(stateR, "el");
  SummStats(s,
            length(VECTOR_ELT(elR, 0)),
            (Vertex*) INTEGER(VECTOR_ELT(elR, 0)),
            (Vertex*) INTEGER(VECTOR_ELT(elR, 1)),
            REAL(stats));

  PutRNGstate();
  UNPROTECT(1);
  return stats;
}


/****************
 void SummStats Computes summary statistics for a network. Must be
 passed an empty network and passed an empty network
*****************/
void SummStats(ErgmState *s, Edge n_edges, Vertex *tails, Vertex *heads, double *stats){
  Network *nwp = s->nwp;
  Model *m = s->m;

  DetShuffleEdges(tails,heads,n_edges); /* Shuffle edgelist. */
  
  Edge ntoggles = n_edges; // So that we can use the macros

  /* Calculate statistics for terms that don't have c_functions or s_functions.  */
  EXEC_THROUGH_TERMS_INTO(m, stats, {
      if(mtp->s_func==NULL && mtp->c_func==NULL && mtp->d_func){
	(*(mtp->d_func))(ntoggles, tails, heads,
			 mtp, nwp);  /* Call d_??? function */
        addonto(dstats, mtp->dstats, N_CHANGE_STATS);
      }
    });

  /* Calculate statistics for terms that have c_functions but not s_functions.  */
  for(Edge e=0; e<n_edges; e++){
    Vertex t=TAIL(e), h=HEAD(e);
    Rboolean edgeflag = IS_OUTEDGE(t,h);

    ergm_PARALLEL_FOR_LIMIT(m->n_terms)        
    EXEC_THROUGH_TERMS_INTO(m, stats, {
	if(mtp->s_func==NULL && mtp->c_func){
	  ZERO_ALL_CHANGESTATS();
	  (*(mtp->c_func))(t, h,
			   mtp, nwp, edgeflag);  /* Call c_??? function */
          addonto(dstats, mtp->dstats, N_CHANGE_STATS);
	}
      });
    
    /* Update storage and network */    
    UPDATE_STORAGE_COND(t, h, nwp, m, NULL, edgeflag, mtp->s_func==NULL && mtp->d_func==NULL);
    TOGGLE_KNOWN(t, h, edgeflag);
  }
  
  /* Calculate statistics for terms have s_functions  */
  EXEC_THROUGH_TERMS_INTO(m, stats, {
      if(mtp->s_func){
	ZERO_ALL_CHANGESTATS();
	(*(mtp->s_func))(mtp, nwp);  /* Call d_??? function */
	for(unsigned int k=0; k<N_CHANGE_STATS; k++){
	  dstats[k] = mtp->dstats[k]; // Overwrite, not accumulate.
	}
      }
    });
}

