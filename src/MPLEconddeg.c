/*  File src/MPLEconddeg.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2014 Statnet Commons
 */
#include "MPLEconddeg.h"

/*****************
 Note on undirected networks:  For j<k, edge {j,k} should be stored
 as (j,k) rather than (k,j).  In other words, only directed networks
 should have (k,j) with k>j.
*****************/

/*****************
 void MPLEconddeg_wrapper

 Wrapper for a call from R.
*****************/
void MPLEconddeg_wrapper(int *dnumnets, int *nedges,
                         int *tails, int *heads,
                   int *dn, int *dflag, int *bipartite, 
                   int *nterms, char **funnames,
                   char **sonames, 
                   char **MHproposaltype, char **MHproposalpackage,
                   double *inputs, double *theta0, int *samplesize, 
                   double *sample, int *burnin, int *interval,  
                   int *newnetworktails, 
                   int *newnetworkheads, 
                   int *fVerbose, 
                   int *attribs, int *maxout, int *maxin, int *minout,
                   int *minin, int *condAllDegExact, int *attriblength, 
                   int *maxedges,
                   int *status){
  int directed_flag;
  Vertex n_nodes, nmax, bip;
  /* Edge n_networks; */
  Network nw[1];
  Model *m;
  MHproposal MH;
  
  n_nodes = (Vertex)*dn; /* coerce double *dn to type Vertex */
  /* n_networks = (Edge)*dnumnets; */ /* coerce double *dnedges to type Edge */
  nmax = (Edge)*maxedges; /* coerce double *maxedges to type Edge */
  bip = (Vertex)*bipartite; /* coerce double *bipartite to type Vertex */
  
  GetRNGstate();  /* R function enabling uniform RNG */
  
  directed_flag = *dflag;

  m=ModelInitialize(*funnames, *sonames, &inputs, *nterms);

  /* Form the missing network */
  nw[0]=NetworkInitialize(tails, heads, nedges[0], 
                          n_nodes, directed_flag, bip, 0, 0, NULL);
 
  MH_init(&MH,
	  *MHproposaltype, *MHproposalpackage,
	  inputs,
	  *fVerbose,
	  nw, attribs, maxout, maxin, minout, minin,
	  *condAllDegExact, *attriblength);

  *status = CondDegSampler (&MH,
	      theta0, sample, *samplesize,
	      *burnin, *interval,
	      *fVerbose, nmax, nw, m);

  MH_free(&MH);
        
/* Rprintf("Back! %d %d\n",nw[0].nedges, nmax); */

  /* record new generated network to pass back to R */
  if(*status == MCMC_OK && *maxedges>0 && newnetworktails && newnetworkheads)
    newnetworktails[0]=newnetworkheads[0]=EdgeTree2EdgeList(newnetworktails+1,newnetworkheads+1,nw,nmax-1);
  
  ModelDestroy(m);

  NetworkDestroy(nw);
  PutRNGstate();  /* Disable RNG before returning */
}


/*********************
 MCMCStatus CondDegSampler

 Using the parameters contained in the array theta, obtain the
 network statistics for a sample of size samplesize.  burnin is the
 initial number of Markov chain steps before sampling anything
 and interval is the number of MC steps between successive 
 networks in the sample.  Put all the sampled statistics into
 the networkstatistics array. 
*********************/
MCMCStatus CondDegSampler (MHproposal *MHp,
  double *theta, double *networkstatistics, 
  int samplesize, int burnin, 
  int interval, int fVerbose, int nmax,
  Network *nwp, Model *m) {
//int staken, tottaken, j;
  int i;
  
  /*********************
  networkstatistics are modified in groups of m->n_stats, and they
  reflect the CHANGE in the values of the statistics from the
  original (observed) network.  Thus, when we begin, the initial 
  values of the first group of m->n_stats networkstatistics should 
  all be zero
  *********************/

  /*********************
   Burn in step.  While we're at it, use burnin statistics to 
   prepare covariance matrix for Mahalanobis distance calculations 
   in subsequent calls to M-H
   *********************/
  /*  Catch more edges than we can return */
 //     Rprintf("Sample size %d\n", samplesize);
//      Rprintf("burnin %d interval %d\n", burnin, interval);
//  if(CondDegSample(MHp, theta, networkstatistics, burnin, &staken,
//                   fVerbose, nwp, m)!=MCMC_OK)
//     return MCMC_MH_FAILED;
//     if(nmax!=0 && nwp->nedges >= nmax-1){
//       return MCMC_TOO_MANY_EDGES;
//     }
/*   if (fVerbose){ 
       Rprintf(".");
     } */
  //    Rprintf("From %f theta0\n", 333);
  MHp->logratio = 0;
  
  if (samplesize>1){
//  staken = 0;
//  tottaken = 0;
    
    /* Now sample networks */
    for (i=1; i < samplesize; i++){
      networkstatistics += m->n_stats;
      /* This then adds the change statistics to these values */
      
      (*(MHp->func))(MHp, nwp); /* Call MH function to propose toggles */
    
      /* Calculate change statistics. */
      ChangeStats(MHp->ntoggles, MHp->toggletail, MHp->togglehead, nwp, m);
      
      /* record network statistics for posterity */
      for (unsigned int j = 0; j < m->n_stats; j++){
        networkstatistics[j] += m->workspace[j];
      }
//    Rprintf("n_stats %d ns %f %f\n", m->n_stats,  networkstatistics[0],networkstatistics[1]);

#ifdef Win32
      if( ((100*i) % samplesize)==0 && samplesize > 500){
	R_FlushConsole();
    	R_ProcessEvents();
      }
#endif
  //    Rprintf("Sample size %d\n", samplesize);
  //    Rprintf("Sampled %d from CondDegSample\n", i);
      /* if (fVerbose){
        if( ((3*i) % samplesize)==0 && samplesize > 500){
        Rprintf("Sampled %d from CondDegSample\n", i);}
      }
      if( ((3*i) % samplesize)==0 && tottaken == ptottaken){
        ptottaken = tottaken; 
        Rprintf("Warning:  CondDegSample algorithm has accepted only "
        "%d steps out of a possible %d\n",  ptottaken-tottaken, i); 
      } 
      if (fVerbose && (i % dotinterval)==0) { 
        Rprintf(".");  
      } */
    }
  }
  return MCMC_OK;
}
