/*  File src/godfather.c in package tergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2008-2017 Statnet Commons
 */
#include "wtMCMC.h"
#include "wtchangestat.h"

WtMCMCStatus WtGodfather(Edge n_changes, Vertex *tails, Vertex *heads, double *weights,
	       WtNetwork *nwp, WtModel *m, double *stats){

  stats+=m->n_stats;
  
  for (unsigned int termi=0; termi < m->n_terms; termi++)
    m->termarray[termi].dstats = m->workspace;
  
  /* Doing this one change at a time saves a lot of changes... */
  for(Edge e=0; e<n_changes; e++){
    WtModelTerm *mtp = m->termarray;
    double *statspos=stats;
    int tail = tails[e], head = heads[e];
    double weight = weights[e];

    if(tail==0){
      stats+=m->n_stats;
      continue;
    }
    
    for (unsigned int termi=0; termi < m->n_terms; termi++, mtp++){
      (*(mtp->d_func))(1, &tail, &head, &weight,
		       mtp, nwp);  /* Call d_??? function */
      for (unsigned int i=0; i < mtp->nstats; i++,statspos++)
	*statspos += mtp->dstats[i];
    }
    
    SETWT(tail,head,weight);
  }

  return WtMCMC_OK;
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
void WtGodfather_wrapper(int *n_edges, int *tails, int *heads, double *weights,
			 int *n_nodes, int *dflag, int *bipartite, 
			 int *nterms, char **funnames, char **sonames, double *inputs,
			 int *total_changes, int *changetails, int *changeheads, double *changeweights,
			 double *changestats, 
			 int *maxedges,
			 int *newnetworktails, 
			 int *newnetworkheads, 
			 double *newnetworkweights, 
			 int *fVerbose, 
			 int *status){
  Vertex nmax;
  /* Edge n_networks; */
  WtNetwork nw[1];
  WtModel *m;
  
  /* n_networks = (Edge)*dnumnets;  */
  nmax = (Edge)abs(*maxedges);

  m=WtModelInitialize(*funnames, *sonames, &inputs, *nterms);

  /* Form the network */
  nw[0]=WtNetworkInitialize(tails, heads, weights, n_edges[0], 
			    *n_nodes, *dflag, *bipartite, 0, 0, NULL);
  
  *status = WtGodfather(abs(*total_changes), changetails, changeheads, changeweights,
			nw, m, changestats);
  
  /* record new generated network to pass back to R */
  if(*status == WtMCMC_OK && *maxedges>0 && newnetworktails && newnetworkheads && newnetworkweights)
    newnetworktails[0]=newnetworkheads[0]=WtEdgeTree2EdgeList(newnetworktails+1,newnetworkheads+1,newnetworkweights+1,nw,nmax-1);
  
  WtModelDestroy(m);
  WtNetworkDestroy(nw);
}
