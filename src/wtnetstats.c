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
/*****************
 void network_stats_wrapper

 Wrapper for a call from R.  Return the change in the statistics when
 we go from an empty graph to the observed graph.  If the empty graph
 has true global values equal to zero for all statistics, then this
 change gives the true global values for the observed graph.
*****************/
void wt_network_stats_wrapper(int *tails, int *heads, double *weights, int *timings, int *time, int *lasttoggle, int *dnedges,
			   int *dn, int *dflag,  int *bipartite,
			   int *nterms, char **funnames,
			   char **sonames, double *inputs,  double *stats)
{
  int directed_flag;
  Vertex n_nodes;
  Edge n_edges;
  WtNetwork *nwp;
  WtModel *m;
  Vertex bip;

/*	     Rprintf("prestart with setup\n"); */
  n_nodes = (Vertex)*dn; 
  n_edges = (Edge)*dnedges;     
  directed_flag = *dflag;
  bip = (Vertex)*bipartite;
  
  if(*lasttoggle == 0) lasttoggle = NULL;

  GetRNGstate();  /* R function enabling uniform RNG */

  m=WtModelInitialize(*funnames, *sonames, &inputs, *nterms);
  nwp=WtNetworkInitialize(NULL, NULL, NULL, 0,
			    n_nodes, directed_flag, bip, *timings?1:0, *timings?*time:0, *timings?lasttoggle:NULL);

  /* Compute the change statistics and copy them to stats for return
     to R.  Note that stats already has the statistics of an empty
     network, so d_??? statistics will add on to them, while s_???
     statistics will simply overwrite them.*/
  WtSummStats(n_edges, (Vertex*)tails, (Vertex*)heads, weights, nwp, m,stats);
  
  WtModelDestroy(nwp, m);
  WtNetworkDestroy(nwp);
  PutRNGstate();
}


/****************
 void SummStats Computes summary statistics for a network. Must be
 passed an empty network (and a possible discordance network) and 
 passed an empty network
*****************/
void WtSummStats(Edge n_edges, Vertex *tails, Vertex *heads, double *weights,
WtNetwork *nwp, WtModel *m, double *stats){

  WtDetShuffleEdges(tails,heads,weights,n_edges); /* Shuffle edgelist. */
  
  Edge ntoggles = n_edges; // So that we can use the macros

  /* Initialize storage for terms that don't have s_functions.  */
  WtEXEC_THROUGH_TERMS_INREVERSE(m, {
      IFDEBUG_BACKUP_DSTATS;
      if(mtp->s_func==NULL && mtp->i_func)
	(*(mtp->i_func))(mtp, nwp);  /* Call i_??? function */
      else if(mtp->s_func==NULL && mtp->u_func) /* No initializer but an updater -> uses a 1-function implementation. */
	(*(mtp->u_func))(0, 0, 0, mtp, nwp);  /* Call u_??? function */
      IFDEBUG_RESTORE_DSTATS;
    });
    
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
    GETNEWTOGGLEINFO();
    
    ergm_PARALLEL_FOR_LIMIT(m->n_terms)
    WtEXEC_THROUGH_TERMS_INTO(m, stats, {
	if(mtp->s_func==NULL && mtp->c_func){
	  ZERO_ALL_CHANGESTATS();
	  (*(mtp->c_func))(TAIL, HEAD, NEWWT,
			   mtp, nwp);  /* Call c_??? function */
	  
	  for(unsigned int k=0; k<N_CHANGE_STATS; k++){
	    dstats[k] += mtp->dstats[k];
	  }
	}
      });
    
    /* Update storage and network */    
    WtUPDATE_STORAGE_COND(TAIL, HEAD, NEWWT, nwp, m, NULL, mtp->s_func==NULL && mtp->d_func==NULL);
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

