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

 Wrapper for a call from R.  Return the change in the statistics when
 we go from an empty graph to the observed graph.  If the empty graph
 has true global values equal to zero for all statistics, then this
 change gives the true global values for the observed graph.
*****************/

SEXP wt_network_stats_wrapper(ARGS_WTNWSETTINGS,
                              ARGS_WTMODEL,
                              ARGS_WTINPUTS,
                              ARGS_WTNWSTATE,
                              ARGS_WTLASTTOGGLE,
                              // Summary settings
                              SEXP emptynwstats){
  GetRNGstate();  /* R function enabling uniform RNG */
  WtErgmState *s = WtErgmStateInit(YES_WTNWSETTINGS,
                                   YES_WTMODEL_NOINIT_S,
                                   NO_WTMHPROPOSAL,
                                   YES_WTINPUTS,
                                   NO_WTNWSTATE,
                                   YES_WTLASTTOGGLE);

  WtModel *m = s->m;

  SEXP stats = PROTECT(allocVector(REALSXP, m->n_stats));

  if(length(emptynwstats)>0) memcpy(REAL(stats), REAL(emptynwstats), m->n_stats*sizeof(double));
  else memset(REAL(stats), 0, m->n_stats*sizeof(double));

  /* Compute the change statistics and copy them to stats for return
     to R.  Note that stats already has the statistics of an empty
     network, so d_??? statistics will add on to them, while s_???
     statistics will simply overwrite them.*/
  WtSummStats(s, asInteger(nedges), (Vertex*)INTEGER(tails), (Vertex*)INTEGER(heads), REAL(weights), REAL(stats));
  
  PutRNGstate();
  UNPROTECT(1);
  return stats;
}


/****************
 void SummStats Computes summary statistics for a network. Must be
 passed an empty network and passed an empty network
*****************/
void WtSummStats(WtErgmState *s, Edge n_edges, Vertex *tails, Vertex *heads, double *weights, double *stats){
  WtNetwork *nwp = s->nwp;
  WtModel *m = s->m;

  WtDetShuffleEdges(tails,heads,weights,n_edges); /* Shuffle edgelist. */
  
  Edge ntoggles = n_edges; // So that we can use the macros

  /* Calculate statistics for terms that don't have c_functions or s_functions.  */
  WtEXEC_THROUGH_TERMS_INTO(m, stats, {
      if(mtp->s_func==NULL && mtp->c_func==NULL && mtp->d_func){
	(*(mtp->d_func))(ntoggles, tails, heads, weights,
			 mtp, nwp);  /* Call d_??? function */
	for(unsigned int k=0; k<N_CHANGE_STATS; k++){
	  dstats[k] += mtp->dstats[k];
	}
      }
    });

  /* Calculate statistics for terms that have c_functions but not s_functions.  */
  FOR_EACH_TOGGLE{
    GETTOGGLEINFO();
    
    ergm_PARALLEL_FOR_LIMIT(m->n_terms)
    WtEXEC_THROUGH_TERMS_INTO(m, stats, {
	if(mtp->s_func==NULL && mtp->c_func){
	  ZERO_ALL_CHANGESTATS();
	  (*(mtp->c_func))(TAIL, HEAD, NEWWT,
			   mtp, nwp, OLDWT);  /* Call c_??? function */
	  
	  for(unsigned int k=0; k<N_CHANGE_STATS; k++){
	    dstats[k] += mtp->dstats[k];
	  }
	}
      });
    
    /* Update storage and network */    
    WtUPDATE_STORAGE_COND(TAIL, HEAD, NEWWT, nwp, m, NULL, OLDWT, mtp->s_func==NULL && mtp->d_func==NULL);
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

