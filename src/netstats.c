/*  File src/netstats.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2014 Statnet Commons
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
                          n_nodes, directed_flag, bip, *timings?1:0, *timings?*time:0, *timings?lasttoggle:NULL);

  /* Compute the change statistics and copy them to stats for return
     to R.  Note that stats already has the statistics of an empty
     network, so d_??? statistics will add on to them, while s_???
     statistics will simply overwrite them. */
  SummStats(n_edges, tails, heads, nw, m, stats);
  
  ModelDestroy(m);
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
  
  for (unsigned int termi=0; termi < m->n_terms; termi++)
    m->termarray[termi].dstats = m->workspace;
  
  /* Doing this one toggle at a time saves a lot of toggles... */
  for(Edge e=0; e<n_edges; e++){
    ModelTerm *mtp = m->termarray;
    double *statspos=stats;
    
    for (unsigned int termi=0; termi < m->n_terms; termi++, mtp++){
      if(!mtp->s_func){
        (*(mtp->d_func))(1, tails+e, heads+e, 
        mtp, nwp);  /* Call d_??? function */
        for (unsigned int i=0; i < mtp->nstats; i++,statspos++)
          *statspos += mtp->dstats[i];
      }else statspos += mtp->nstats;
    }
    
    ToggleEdge(tails[e],heads[e],nwp);
  }
  
  ModelTerm *mtp = m->termarray;
  double *dstats = m->workspace;
  double *statspos=stats;
  for (unsigned int termi=0; termi < m->n_terms; termi++, dstats+=mtp->nstats, mtp++ ){
    if(mtp->s_func){
      (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */
      for (unsigned int i=0; i < mtp->nstats; i++,statspos++)
        *statspos = mtp->dstats[i];
    }else statspos += mtp->nstats;
  }
  
  PutRNGstate();
}

