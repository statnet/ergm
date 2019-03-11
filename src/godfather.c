/*  File src/godfather.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2019 Statnet Commons
 */
#include "MCMC.h"
#include "ergm_model.h"
#include "ergm_changestat.h"

MCMCStatus Godfather(Edge n_changes, Vertex *tails, Vertex *heads, int *weights,
	       Network *nwp, Model *m, double *stats){

  stats+=m->n_stats;

  /* Doing this one change at a time saves a lot of changes... */
  for(Edge e=0; e<n_changes; e++){
    Vertex t=TAIL(e), h=HEAD(e); 

    if(t==0){
      stats+=m->n_stats;
      continue;
    }
    
    if(weights){
      if(IS_OUTEDGE(t,h)==weights[e])
	continue;
    }

    EXEC_THROUGH_TERMS_INTO(m, stats, {
	if(mtp->c_func){
	  ZERO_ALL_CHANGESTATS();
	  (*(mtp->c_func))(t, h,
			   mtp, nwp);  /* Call c_??? function */
	}else if(mtp->d_func){
	  (*(mtp->d_func))(1, &t, &h,
			   mtp, nwp);  /* Call d_??? function */
	}
	for(unsigned int k=0; k<N_CHANGE_STATS; k++){
	    dstats[k] += mtp->dstats[k];
	}
      });


    /* Update storage and network */    
    UPDATE_STORAGE(t, h, nwp, m, NULL);
    TOGGLE(t,h);
  } 

  return MCMC_OK;
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
void Godfather_wrapper(int *n_edges, int *tails, int *heads,
		       int *n_nodes, int *dflag, int *bipartite, 
		       int *nterms, char **funnames, char **sonames, double *inputs,
		       int *total_changes, int *changetails, int *changeheads, int *changeweights,
		       double *changestats, 
		       int *maxedges,
		       int *newnetworktails, 
		       int *newnetworkheads, 
		       int *fVerbose, 
		       int *status){
  Vertex nmax;
  /* Edge n_networks; */
  Network *nwp;
  Model *m;
  
  nmax = (Edge)abs(*maxedges);

  GetRNGstate();  /* R function enabling uniform RNG */

  m=ModelInitialize(*funnames, *sonames, &inputs, *nterms);

  /* Form the network */
  nwp=NetworkInitialize((Vertex*)tails, (Vertex*)heads, n_edges[0], 
                          *n_nodes, *dflag, *bipartite, 0, 0, NULL);
  
  /* Trigger initial storage update */
  InitStats(nwp, m);

  *status = Godfather(abs(*total_changes), (Vertex*)changetails, (Vertex*)changeheads, *total_changes<0? NULL : changeweights,
		      nwp, m, changestats);
  
  /* record new generated network to pass back to R */
  if(*status == MCMC_OK && *maxedges>0 && newnetworktails && newnetworkheads)
    newnetworktails[0]=newnetworkheads[0]=EdgeTree2EdgeList((Vertex*)newnetworktails+1,(Vertex*)newnetworkheads+1,nwp,nmax-1);
  
  ModelDestroy(nwp, m);
  NetworkDestroy(nwp);
  PutRNGstate();
}
