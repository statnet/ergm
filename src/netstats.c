/*  File src/netstats.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2013 Statnet Commons
 */
#include "netstats.h"
/*****************
 void network_stats_wrapper

 Wrapper for a call from R.  Return the change in the statistics when
 we go from an empty graph to the observed graph.  If the empty graph
 has true global values equal to zero for all statistics, then this
 change gives the true global values for the observed graph.
*****************/

/* *** don't forget tail-> head, so this fucntion now accepts tails before heads */

void network_stats_wrapper(int *tails, int *heads, int *timings, int *time, int *lasttoggle, int *dnedges, 
			   int *dn, int *dflag,  int *bipartite,
			   int *nterms, char **funnames,
			   char **sonames, double *inputs,  double *stats)
{
  int directed_flag;
  Vertex n_nodes;
  Edge n_edges;
  Network nw[2];
  Model *m;
  Vertex bip;

/*	     Rprintf("prestart with setup\n"); */
  n_nodes = (Vertex)*dn; 
  n_edges = (Edge)*dnedges;     
  directed_flag = *dflag;
  bip = (Vertex)*bipartite;
  
  if(*lasttoggle == 0) lasttoggle = NULL;

  m=ModelInitialize(*funnames, *sonames, &inputs, *nterms);
  nw[0]=NetworkInitialize(NULL, NULL, 0,
                          n_nodes, directed_flag, bip, *timings?1:0, *timings?*time:0, *timings?lasttoggle:NULL, m->n_aux);

  /* Compute the change statistics and copy them to stats for return
     to R.  Note that stats already has the statistics of an empty
     network, so d_??? statistics will add on to them, while s_???
     statistics will simply overwrite them. */
  SummStats(n_edges, tails, heads, nw, m, stats);
  
  ModelDestroy(m, nw);
  NetworkDestroy(nw);
}


/****************
 void SummStats Computes summary statistics for a network. Must be
 passed an empty network (and a possible discordance network) and 
 passed an empty network
*****************/

/* *** don't forget tail-> head, so this fucntion now accepts tails before heads */

void SummStats(Edge n_edges, Vertex *tails, Vertex *heads,
Network *nwp, Model *m, double *stats){
  
  GetRNGstate();  /* R function enabling uniform RNG */
  
  ShuffleEdges(tails,heads,n_edges); /* Shuffle edgelist. */
  
  Edge ntoggles = n_edges; // So that we can use the macros

  /* Initialize storage for terms that don't have s_functions.  */
  EXEC_THROUGH_TERMS({
#ifdef DEBUG
      double *dstats = mtp->dstats;
      mtp->dstats = NULL; // Trigger segfault if i_func tries to write to change statistics.
#endif
      if(mtp->s_func==NULL && mtp->i_func)
	(*(mtp->i_func))(mtp, nwp);  /* Call i_??? function */
      else if(mtp->s_func==NULL && mtp->u_func) /* No initializer but an updater -> uses a 1-function implementation. */
	(*(mtp->u_func))(0, 0, mtp, nwp);  /* Call u_??? function */
#ifdef DEBUG
      mtp->dstats = dstats;
#endif
    });
    
  /* Calculate statistics for terms that don't have c_functions or s_functions.  */
  EXEC_THROUGH_TERMS_INTO(stats, {
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
    
    EXEC_THROUGH_TERMS_INTO(stats, {
	if(mtp->s_func==NULL && mtp->c_func){
	  (*(mtp->c_func))(t, h,
			   mtp, nwp);  /* Call c_??? function */
	  
	  for(unsigned int k=0; k<N_CHANGE_STATS; k++){
	    dstats[k] += mtp->dstats[k];
	  }
	}
      });
    
    /* Update storage and network */    
    UPDATE_C_STORAGE(t, h, m, nwp);
    TOGGLE(t, h);
  }
  
  /* Calculate statistics for terms have s_functions  */
  EXEC_THROUGH_TERMS_INTO(stats, {
      if(mtp->s_func){
	(*(mtp->s_func))(mtp, nwp);  /* Call d_??? function */
	for(unsigned int k=0; k<N_CHANGE_STATS; k++){
	  dstats[k] = mtp->dstats[k]; // Overwrite, not accumulate.
	}
      }
    });

  PutRNGstate();
}

