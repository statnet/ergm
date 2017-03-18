/*  File src/SAN.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2017 Statnet Commons
 */
#include "SAN.h"


/*****************
 Note on undirected networks:  For j<k, edge {j,k} should be stored
 as (j,k) rather than (k,j).  In other words, only directed networks
 should have (k,j) with k>j.
*****************/

/*****************
 void SAN_wrapper

 Wrapper for a call from R.
*****************/

/* *** don't forget tail-> head, so this function now accepts tails before heads */

void SAN_wrapper ( int *dnumnets, int *nedges,
		   int *tails, int *heads,
                   int *dn, int *dflag, int *bipartite, 
                   int *nterms, char **funnames,
                   char **sonames, 
                   char **MHproposaltype, char **MHproposalpackage,
                   double *inputs, double *theta0, double *tau, 
		   int *samplesize, 
                   double *sample, int *burnin, int *interval,  
                   int *newnetworktails, 
                   int *newnetworkheads, 
                   double *invcov, 
                   int *fVerbose, 
                   int *attribs, int *maxout, int *maxin, int *minout,
                   int *minin, int *condAllDegExact, int *attriblength, 
                   int *maxedges,
		   int *status){
  int directed_flag;
  Vertex n_nodes, nmax, bip;
  Network nw[1];
  Model *m;
  MHproposal MH;


  /* please don't forget:   tail -> head   */

  
  n_nodes = (Vertex)*dn; 
  nmax = (Edge)abs(*maxedges);
  bip = (Vertex)*bipartite; 
  
  GetRNGstate();  /* R function enabling uniform RNG */
  
  directed_flag = *dflag;

  m=ModelInitialize(*funnames, *sonames, &inputs, *nterms);

  /* Form the network */
  nw[0]=NetworkInitialize(tails, heads, nedges[0], 
                          n_nodes, directed_flag, bip, 0, 0, NULL);
  
  MH_init(&MH,
	  *MHproposaltype, *MHproposalpackage,
	  inputs,
	  *fVerbose,
	  nw, attribs, maxout, maxin, minout, minin,
	  *condAllDegExact, *attriblength);

  *status = SANSample (&MH,
		       theta0, invcov, tau, sample, *samplesize,
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
 MCMCStatus SANSample

 Using the parameters contained in the array theta, obtain the
 network statistics for a sample of size samplesize.  burnin is the
 initial number of Markov chain steps before sampling anything
 and interval is the number of MC steps between successive 
 networks in the sample.  Put all the sampled statistics into
 the networkstatistics array. 
*********************/
MCMCStatus SANSample (MHproposal *MHp,
  double *theta, double *invcov, double *tau, double *networkstatistics, 
  int samplesize, int burnin, 
  int interval, int fVerbose, int nmax,
  Network *nwp, Model *m) {
  int staken, tottaken, ptottaken;
    
  /*********************
  networkstatistics are modified in groups of m->n_stats, and they
  reflect the CHANGE in the values of the statistics from the
  original (observed) network.  Thus, when we begin, the initial 
  values of the first group of m->n_stats networkstatistics should 
  all be zero
  *********************/
/*for (j=0; j < m->n_stats; j++) */
/*  networkstatistics[j] = 0.0; */
/* Rprintf("\n"); */
/* for (j=0; j < m->n_stats; j++){ */
/*   Rprintf("j %d %f\n",j,networkstatistics[j]); */
/* } */
/* Rprintf("\n"); */

  /*********************
   Burn in step.  While we're at it, use burnin statistics to 
   prepare covariance matrix for Mahalanobis distance calculations 
   in subsequent calls to M-H
   *********************/
  /*  Catch more edges than we can return */
  if(SANMetropolisHastings(MHp, theta, invcov, tau, networkstatistics, burnin, &staken,
			   fVerbose, nwp, m)!=MCMC_OK)
    return MCMC_MH_FAILED;
  if(nmax!=0 && nwp->nedges >= nmax-1){
    return MCMC_TOO_MANY_EDGES;
  }
  
  if (fVerbose){
    Rprintf("Returned from SAN Metropolis-Hastings burnin\n");
  }
  
  if (samplesize>1){
    staken = 0;
    tottaken = 0;
    ptottaken = 0;
    
    /* Now sample networks */
    for (unsigned int i=1; i < samplesize; i++){
      /* Set current vector of stats equal to previous vector */
      for (unsigned int j=0; j<m->n_stats; j++){
        networkstatistics[j+m->n_stats] = networkstatistics[j];
      }
      networkstatistics += m->n_stats;
      /* This then adds the change statistics to these values */
      
      if(SANMetropolisHastings(MHp, theta, invcov, tau, networkstatistics, 
		             interval, &staken,
			       fVerbose, nwp, m)!=MCMC_OK)
	return MCMC_MH_FAILED;
      if(nmax!=0 && nwp->nedges >= nmax-1){
	return MCMC_TOO_MANY_EDGES;
      }
      tottaken += staken;
      if (fVerbose){
        if( ((3*i) % samplesize)==0 && samplesize > 500){
        Rprintf("Sampled %d from SAN Metropolis-Hastings\n", i);}
      }
      
      if( ((3*i) % samplesize)==0 && tottaken == ptottaken){
        ptottaken = tottaken; 
        Rprintf("Warning:  SAN Metropolis-Hastings algorithm has accepted only "
        "%d steps out of a possible %d\n",  ptottaken-tottaken, i); 
      }

#ifdef Win32
      if( ((100*i) % samplesize)==0 && samplesize > 500){
	R_FlushConsole();
    	R_ProcessEvents();
      }
#endif
    }
    /*********************
    Below is an extremely crude device for letting the user know
    when the chain doesn't accept many of the proposed steps.
    *********************/
    if (fVerbose){
	  if (samplesize > 0 && interval > LONG_MAX / samplesize) {
		// overflow
		Rprintf("SAN Metropolis-Hastings accepted %7.3f%% of %d proposed steps.\n",
	      tottaken*100.0/(1.0*interval*samplesize), interval, samplesize); 
	  } else {
	    Rprintf("SAN Metropolis-Hastings accepted %7.3f%% of %d proposed steps.\n",
	      tottaken*100.0/(1.0*interval*samplesize), interval*samplesize); 
	  }
    }
  }else{
    if (fVerbose){
      Rprintf("SAN Metropolis-Hastings accepted %7.3f%% of %d proposed steps.\n",
	      staken*100.0/(1.0*burnin), burnin); 
    }
  }
  return MCMC_OK;
}

/*********************
MCMCStatus SANMetropolisHastings

 In this function, theta is a m->n_stats-vector just as in SANSample,
 but now networkstatistics is merely another m->n_stats-vector because
 this function merely iterates nsteps times through the Markov
 chain, keeping track of the cumulative change statistics along
 the way, then returns, leaving the updated change statistics in
 the networkstatistics vector.  In other words, this function 
 essentially generates a sample of size one
*********************/
MCMCStatus SANMetropolisHastings (MHproposal *MHp,
			    double *theta, double *invcov, 
			    double *tau, double *networkstatistics,
			    int nsteps, int *staken,
			    int fVerbose,
			    Network *nwp,
			    Model *m) {
  unsigned int taken=0, unsuccessful=0;
  double *deltainvsig, *delta;
  deltainvsig = (double *)malloc( m->n_stats * sizeof(double));
  delta = (double *)malloc( m->n_stats * sizeof(double));
  
/*  if (fVerbose)
    Rprintf("Now proposing %d MH steps... ", nsteps); */
  for(unsigned int step=0; step < nsteps; step++) {
    MHp->logratio = 0;
    (*(MHp->func))(MHp, nwp); /* Call MH function to propose toggles */

      if(MHp->toggletail[0]==MH_FAILED){
	switch(MHp->togglehead[0]){
	case MH_UNRECOVERABLE:
	  error("Something very bad happened during proposal. Memory has not been deallocated, so restart R soon.");
	  
	case MH_IMPOSSIBLE:
	  Rprintf("MH Proposal function encountered a configuration from which no toggle(s) can be proposed.\n");
	  return MCMC_MH_FAILED;
	  
	case MH_UNSUCCESSFUL:
	  warning("MH Proposal function failed to find a valid proposal.");
	  unsuccessful++;
	  if(unsuccessful>taken*MH_QUIT_UNSUCCESSFUL){
	    Rprintf("Too many MH Proposal function failures.\n");
	    return MCMC_MH_FAILED;
	  }
	case MH_CONSTRAINT:
	  continue;
      }
    }
    
    if(fVerbose>=5){
      Rprintf("Proposal: ");
      for(unsigned int i=0; i<MHp->ntoggles; i++)
	Rprintf(" (%d, %d)", MHp->toggletail[i], MHp->togglehead[i]);
      Rprintf("\n");
    }

    /* Calculate change statistics,
       remembering that tail -> head */
    ChangeStats(MHp->ntoggles, MHp->toggletail, MHp->togglehead, nwp, m);

    if(fVerbose>=5){
      Rprintf("Changes: (");
      for(unsigned int i=0; i<m->n_stats; i++)
	Rprintf(" %f ", m->workspace[i]);
      Rprintf(")\n");
    }
    
    /* Calculate inner product */
    double ip=0, dif=0;
    for (unsigned int i=0; i<m->n_stats; i++){
     delta[i]=0.0;
     deltainvsig[i]=0.0;
     for (unsigned int j=0; j<m->n_stats; j++){
      delta[i]+=networkstatistics[j]*invcov[i+(m->n_stats)*j];
      deltainvsig[i]+=(m->workspace[j])*invcov[i+(m->n_stats)*j];
     }
     ip+=deltainvsig[i]*((m->workspace[i])+2.0*networkstatistics[i]);
     dif+=delta[i]*networkstatistics[i];
    }
    if(fVerbose>=5){
      Rprintf("log acceptance probability: %f\n", ip);
    }
    
    /* if we accept the proposed network */
    if (ip <= 0.0) { 
//  if (ip <= 0.0 || (ip/dif) < 0.001) { 
//  if (div > 0.0 && (ip < 0.0 || unif_rand() < 0.01)) { 
// if (ip <= 0.0 || (ip/dif) < (nsteps-step)*0.001*tau[0]/(1.0*nsteps)) { 
//  if (ip > exp(theta[0])*(m->n_stats)*unif_rand()/(1.0+exp(theta[0])) { 
//  if (ip > tau[0]*(m->n_stats)*unif_rand()) { 
      if(fVerbose>=5){
	Rprintf("Accepted.\n");
      }

      /* Make proposed toggles (updating timestamps--i.e., for real this time) */
      for(unsigned int i=0; i < MHp->ntoggles; i++){
	ToggleEdge(MHp->toggletail[i], MHp->togglehead[i], nwp);
	
	if(MHp->discord)
	  for(Network **nwd=MHp->discord; *nwd!=NULL; nwd++){
	    ToggleEdge(MHp->toggletail[i],  MHp->togglehead[i], *nwd);
	  }
      }
      /* record network statistics for posterity */
      for (unsigned int i = 0; i < m->n_stats; i++){
	networkstatistics[i] += m->workspace[i];
      }
      taken++;
    }else{
      if(fVerbose>=5){
	Rprintf("Rejected.\n");
      }
    }
  }

  free(deltainvsig);
  free(delta);

  *staken = taken;
  return MCMC_OK;
}
