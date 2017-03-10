/*  File src/wtnetstats.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2013 Statnet Commons
 */
#include "wtnetstats.h"
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
  WtNetwork nw[2];
  WtModel *m;
  Vertex bip;

/*	     Rprintf("prestart with setup\n"); */
  n_nodes = (Vertex)*dn; 
  n_edges = (Edge)*dnedges;     
  directed_flag = *dflag;
  bip = (Vertex)*bipartite;
  
  if(*lasttoggle == 0) lasttoggle = NULL;

  m=WtModelInitialize(*funnames, *sonames, &inputs, *nterms);
  nw[0]=WtNetworkInitialize(NULL, NULL, NULL, 0,
			    n_nodes, directed_flag, bip, *timings?1:0, *timings?*time:0, *timings?lasttoggle:NULL, m->n_aux);

  /* Compute the change statistics and copy them to stats for return
     to R.  Note that stats already has the statistics of an empty
     network, so d_??? statistics will add on to them, while s_???
     statistics will simply overwrite them.*/
  WtSummStats(n_edges, tails, heads, weights, nw, m,stats);
  
  WtModelDestroy(m, nw);
  WtNetworkDestroy(nw);
}


/****************
 void SummStats Computes summary statistics for a network. Must be
 passed an empty network (and a possible discordance network) and 
 passed an empty network
*****************/
void WtSummStats(Edge n_edges, Vertex *tails, Vertex *heads, double *weights,
WtNetwork *nwp, WtModel *m, double *stats){

  GetRNGstate();  /* R function enabling uniform RNG */
  
  WtShuffleEdges(tails,heads,weights,n_edges); /* Shuffle edgelist. */
  
  memset(m->workspace, 0, m->n_stats*sizeof(double)); /* Zero all change stats. */ 

  Edge ntoggles = n_edges; // So that we can use the macros

  /* Initialize storage for terms that don't have s_functions.  */
  WtModelTerm *mtp = m->termarray;
  for (unsigned int i=0; i < m->n_terms; i++){
    double *dstats = mtp->dstats;
    mtp->dstats = NULL; // Trigger segfault if i_func tries to write to change statistics.
    if(mtp->s_func==NULL && mtp->i_func)
      (*(mtp->i_func))(mtp, nwp);  /* Call i_??? function */
    else if(mtp->s_func==NULL && mtp->u_func) /* No initializer but an updater -> uses a 1-function implementation. */
      (*(mtp->u_func))(0, 0, 0, mtp, nwp);  /* Call u_??? function */
    mtp->dstats = dstats;
    mtp++;
  }
  
  /* Calculate statistics for terms that don't have c_functions or s_functions.  */
  mtp = m->termarray;
  double *dstats = m->workspace;
  for (unsigned int i=0; i < m->n_terms; i++){
    mtp->dstats = dstats; /* Stuck the change statistic here.*/
    if(mtp->s_func==NULL && mtp->c_func==NULL && mtp->d_func)
      (*(mtp->d_func))(ntoggles, tails, heads, weights,
		       mtp, nwp);  /* Call d_??? function */
    dstats += (mtp++)->nstats;
  }

  
  /* Make a pass through terms with c_functions. */
  FOR_EACH_TOGGLE{
    GETNEWTOGGLEINFO();
    
    WtModelTerm *mtp = m->termarray;
    double *dstats = m->workspace;

    /* Calculate change statistics */
    for (unsigned int i=0; i < m->n_terms; i++){
      mtp->dstats = m->dstatarray[i]; /* If only one toggle, just write directly into the workspace array. */
      if(mtp->s_func==NULL && mtp->c_func)
	(*(mtp->c_func))(TAIL, HEAD, NEWWT,
			 mtp, nwp);  /* Call d_??? function */

      for(unsigned int k=0; k<N_CHANGE_STATS; k++){
	dstats[k] += mtp->dstats[k];
      }
      
      // Advance to the next term.
      dstats += mtp->nstats;
      mtp++;
    }
    
    /* Update storage and network */    
    UPDATE_C_STORAGE(TAIL, HEAD, NEWWT, m, nwp);
    SETWT(TAIL, HEAD, NEWWT);
  }
  
  /* Calculate statistics for terms have s_functions  */
  mtp = m->termarray;
  dstats = m->workspace;
  for (unsigned int i=0; i < m->n_terms; i++){
    mtp->dstats = dstats; /* Stuck the change statistic here.*/
    if(mtp->s_func)
      (*(mtp->s_func))(mtp, nwp);  /* Call d_??? function */
    dstats += (mtp++)->nstats;
  }

  memcpy(stats, m->workspace, m->n_stats*sizeof(double));
  
  PutRNGstate();
}

