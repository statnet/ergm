/*
 *  File ergm/src/godfather.c
 *  Part of the statnet package, http://statnet.org
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) in
 *    http://statnet.org/attribution
 *
 *  Copyright 2012 the statnet development team
 */
#include "godfather.h"

/*****************
 void godfather_wrapper

 ...we'll make them an offer (of toggles) they can't refuse.
 This function takes a list of toggles, each with a time stamp,
 then produces a matrix of changestats (with one row for each unique
 time stamp value) that result from performing all the toggles at
 each time step.  For instance, one might use this function to 
 find the changestats that result from starting from an empty network
 and then adding all of the edges to make up an observed network of interest.
*****************/
void godfather_wrapper (int *tails, int *heads, int *dnedges,
			int *dn, int *dflag, int *bipartite, 
			int *nterms, char **funnames,
			char **sonames, 
			int *totalntoggles, int *timestamps, 
			int *toggletails, int *toggleheads,
			int *dstart, int *dend,
			double *inputs, 
			double *changestats, 
			int *newnetworktails, 
			int *newnetworkheads, 
			int *accumulate, 
			int *fVerbose, 
			int *maxedges) {
  int directed_flag;
  Vertex n_nodes, bip;
  Edge j, n_edges, nmax, tnt;
  Network nw;
  Model *m;
  int start=*dstart, end=*dend; 
  
  n_nodes = (Vertex)*dn; /* coerce double *dn to type Vertex */
  n_edges = (Edge)*dnedges; /* coerce double *dnedges to type Vertex */
  nmax = (Edge)*maxedges; /* coerce double *maxedges to type Edge */
  bip = (Vertex)*bipartite; /* coerce double *bipartite to type Vertex */    
  tnt = (Edge) *totalntoggles; /* coerce double *totalntoggles to type Edge */ 
  directed_flag = *dflag;

  m=ModelInitialize(*funnames, *sonames, &inputs, *nterms);
  nw = NetworkInitialize(tails, heads, n_edges,
                         n_nodes, directed_flag, bip, 1, 0, NULL);

  /*  if (*fVerbose) {
    Rprintf("Total m->n_stats is %i.\n",
    m->n_stats);
    Rprintf("maxedges = %ld, totalntoggles = %ld\n",
    nmax, tnt);
    }*/
  
  /*********************
  changestats are modified in groups of m->n_stats, and they
  reflect the CHANGE in the values of the statistics from the
  original (observed) network.  Thus, when we begin, the initial 
  values of the first group of m->n_stats changestats should 
  all be zero
  *********************/
  for (j=0; j < m->n_stats; j++)
    changestats[j] = 0.0;
  
  /* Now start obtaining change statistics */
  unsigned int pos=0;
  for(unsigned int t=start;t<=end;t++){
    changestats += m->n_stats;
    for (j=0; j<m->n_stats; j++){
      changestats[j] = changestats[j-m->n_stats];
    }
    for(;pos<tnt && timestamps[pos]==t;pos++){
      ChangeStats(1, toggletails+pos, toggleheads+pos, &nw, m);
      /* Accumulate change statistics */
      for (j=0; j<m->n_stats; j++){
        changestats[j] += m->workspace[j];
      }
	
      /* Make proposed toggles (for real this time) */
      //Obvious bug in next line, since i is uninitialized and j just finished loop:
      //if (!(*accumulate) || EdgetreeSearch(toggletails[i-j], toggleheads[i-j], nw.outedges) == 0) { 
      if (!(*accumulate)) { 
        ToggleEdgeWithTimestamp(toggletails[pos], toggleheads[pos], &nw);
      }
    }
  }

  if (nmax>0) {
    /* record new generated network to pass back to R */
    newnetworktails[0]=newnetworkheads[0]=EdgeTree2EdgeList(newnetworktails+1,newnetworkheads+1,&nw,nmax);
    if (newnetworktails[0]>=nmax) { 
      Rprintf("Error!  The value of maxedges was not set high enough in ergm.godfather\n");
    }
  }
  /* Clean up and return */
  ModelDestroy(m);
  NetworkDestroy(&nw);
}

