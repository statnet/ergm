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

SEXP network_stats_wrapper(// Network settings
                           SEXP dn, SEXP dflag, SEXP bipartite,
                           // Model settings
                           SEXP nterms, SEXP funnames,
                           SEXP sonames,
                           // Numeric inputs
                           SEXP inputs,
                           // Network state
                           SEXP nedges,
                           SEXP tails, SEXP heads,
                           SEXP time, SEXP lasttoggle,
                           // Summary settings
                           SEXP emptynwstats){
  GetRNGstate();  /* R function enabling uniform RNG */
  Rboolean timings = length(time)>0;

  ErgmState *s = ErgmStateInit(// Network settings
                               asInteger(dn), asInteger(dflag), asInteger(bipartite),
                               // Model settings
                               asInteger(nterms), FIRSTCHAR(funnames), FIRSTCHAR(sonames), TRUE,
                               // Proposal settings
                               NO_MHPROPOSAL,
                               // Numeric inputs
                               REAL(inputs),
                               // Network state
                               NO_NWSTATE,
                               timings, asInteger(time), INTEGER(lasttoggle));

  Model *m = s->m;

  SEXP stats = PROTECT(allocVector(REALSXP, m->n_stats));

  if(length(emptynwstats)>0) memcpy(REAL(stats), REAL(emptynwstats), m->n_stats*sizeof(double));
  else memset(REAL(stats), 0, m->n_stats*sizeof(double));

  /* Compute the change statistics and copy them to stats for return
     to R.  Note that stats already has the statistics of an empty
     network, so d_??? statistics will add on to them, while s_???
     statistics will simply overwrite them. */
  SummStats(s, asInteger(nedges), (Vertex*)INTEGER(tails), (Vertex*)INTEGER(heads), REAL(stats));
  
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
	for(unsigned int k=0; k<N_CHANGE_STATS; k++){
	  dstats[k] += mtp->dstats[k];
	}
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
	  
	  for(unsigned int k=0; k<N_CHANGE_STATS; k++){
	    dstats[k] += mtp->dstats[k];
	  }
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

