/*  File src/wtSAN.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2019 Statnet Commons
 */
#include "wtSAN.h"


/*****************
 Note on undirected networks:  For j<k, edge {j,k} should be stored
 as (j,k) rather than (k,j).  In other words, only directed networks
 should have (k,j) with k>j.
*****************/

/*****************
 void WtSAN_wrapper

 Wrapper for a call from R.
*****************/
void WtSAN_wrapper (int *nedges,
		    int *tails, int *heads, double *weights, 
		    int *dn, int *dflag, int *bipartite, 
		    int *nterms, char **funnames,
		    char **sonames, 
		    char **MHProposaltype, char **MHProposalpackage,
		    double *inputs, double *tau, 
		    double *sample, double *prop_sample,
		    int *samplesize, int *nsteps,
		    int *newnetworktails, 
		    int *newnetworkheads, 
		    double *newnetworkweights,
		    double *invcov, 
		    int *fVerbose, 
		    int *maxedges,
		    int *status,
            double *offsets){
  int directed_flag;
  Vertex n_nodes, nmax, bip;
  WtNetwork *nwp;
  WtModel *m;
  WtMHProposal *MHp;
  
  n_nodes = (Vertex)*dn; 
  nmax = (Edge)abs(*maxedges);
  bip = (Vertex)*bipartite; 
  
  GetRNGstate();  /* R function enabling uniform RNG */
  
  directed_flag = *dflag;

  m=WtModelInitialize(*funnames, *sonames, &inputs, *nterms);

  /* Form the network */
  nwp=WtNetworkInitialize((Vertex*)tails, (Vertex*)heads, weights, nedges[0],
			    n_nodes, directed_flag, bip, 0, 0, NULL);

  MHp=WtMHProposalInitialize(
	    *MHProposaltype, *MHProposalpackage,
	    inputs,
	    *fVerbose,
	    nwp);

  if(MHp)
    *status = WtSANSample(MHp,
			  invcov, tau, sample, prop_sample, *samplesize,
			  *nsteps,
			  *fVerbose, nmax, nwp, m, offsets);
  else *status = WtMCMC_MH_FAILED;

  WtMHProposalDestroy(MHp);

/* Rprintf("Back! %d %d\n",nwp[0].nedges, nmax); */

  /* record new generated network to pass back to R */
  /* *** and don't forget edges are (tail, head) */
  if(*status == WtMCMC_OK && *maxedges>0 && newnetworktails && newnetworkheads)
    newnetworktails[0]=newnetworkheads[0]=newnetworkweights[0]=WtEdgeTree2EdgeList((Vertex*)newnetworktails+1,(Vertex*)newnetworkheads+1,newnetworkweights+1,nwp,nmax-1);

  WtModelDestroy(m);
  WtNetworkDestroy(nwp);
  PutRNGstate();  /* Disable RNG before returning */
}


/*********************
 void WtSANSample

 Using the parameters contained in the array theta, obtain the
 network statistics for a sample of size samplesize.  nsteps is the
 initial number of Markov chain steps before sampling anything
 and interval is the number of MC steps between successive 
 networks in the sample.  Put all the sampled statistics into
 the networkstatistics array. 
*********************/
WtMCMCStatus WtSANSample (WtMHProposal *MHp,
  double *invcov, double *tau, double *networkstatistics, double *prop_networkstatistics,
  int samplesize, int nsteps, 
  int fVerbose, int nmax,
  WtNetwork *nwp, WtModel *m,
  double *offsets) {
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

  unsigned int interval = nsteps / samplesize; // Integer division: rounds down.
  unsigned int burnin = nsteps - (samplesize-1)*interval;
  
  /*********************
   Burn in step.  While we're at it, use nsteps statistics to 
   prepare covariance matrix for Mahalanobis distance calculations 
   in subsequent calls to M-H
   *********************/
  /*  Catch more edges than we can return */
  if(WtSANMetropolisHastings(MHp, invcov, tau, networkstatistics, prop_networkstatistics, burnin, &staken,
			     fVerbose, nwp, m, offsets)!=WtMCMC_OK)
    return WtMCMC_MH_FAILED;
  if(nmax!=0 && EDGECOUNT(nwp) >= nmax-1){
    return WtMCMC_TOO_MANY_EDGES;
  }

  if (samplesize>1){
    staken = 0;
    tottaken = 0;
    ptottaken = 0;
    
    /* Now sample networks */
    for (unsigned int i=1; i < samplesize; i++){
      /* Set current vector of stats equal to previous vector */
      Rboolean found = TRUE;
      for (unsigned int j=0; j<m->n_stats; j++){
        if((networkstatistics[j+m->n_stats] = networkstatistics[j])!=0) found = FALSE;
      }
      if(found){
	if(fVerbose) Rprintf("Exact match found.\n");
	break;
      }

      networkstatistics += m->n_stats;
      prop_networkstatistics += m->n_stats;
      /* This then adds the change statistics to these values */
      
      if(WtSANMetropolisHastings (MHp, invcov, tau, networkstatistics, prop_networkstatistics,
		             interval, &staken, fVerbose, nwp, m, offsets)!=WtMCMC_OK)
	return WtMCMC_MH_FAILED;
      if(nmax!=0 && EDGECOUNT(nwp) >= nmax-1){
	return WtMCMC_TOO_MANY_EDGES;
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
	  Rprintf("SAN Metropolis-Hastings accepted %7.3f%% of %lld proposed steps.\n",
	    tottaken*100.0/(1.0*interval*samplesize), (long long) interval*samplesize); 
    }
  }else{
    if (fVerbose){
      Rprintf("SAN Metropolis-Hastings accepted %7.3f%% of %d proposed steps.\n",
	      staken*100.0/(1.0*nsteps), nsteps); 
    }
  }
  return WtMCMC_OK;
}

/*********************
MCMCStatus WtSANMetropolisHastings

 In this function, theta is a m->n_stats-vector just as in SANSample,
 but now networkstatistics is merely another m->n_stats-vector because
 this function merely iterates nsteps times through the Markov
 chain, keeping track of the cumulative change statistics along
 the way, then returns, leaving the updated change statistics in
 the networkstatistics vector.  In other words, this function 
 essentially generates a sample of size one
*********************/
WtMCMCStatus WtSANMetropolisHastings (WtMHProposal *MHp,
			    double *invcov, 
				  double *tau, double *networkstatistics, double *prop_networkstatistics,
			    int nsteps, int *staken,
			    int fVerbose,
			    WtNetwork *nwp,
			    WtModel *m,
                double *offsets) {
  unsigned int taken=0, unsuccessful=0;
  double *deltainvsig;
  deltainvsig = (double *)Calloc(m->n_stats, double);
  
/*  if (fVerbose)
    Rprintf("Now proposing %d WtMH steps... ", nsteps); */
  for(unsigned int step=0; step < nsteps; step++) {
    MHp->logratio = 0;
    (*(MHp->func))(MHp, nwp); /* Call MH function to propose toggles */

    if(MHp->toggletail[0]==MH_FAILED){
      switch(MHp->togglehead[0]){
      case MH_UNRECOVERABLE:
	error("Something very bad happened during proposal. Memory has not been deallocated, so restart R soon.");
	
      case MH_IMPOSSIBLE:
	Rprintf("MH MHProposal function encountered a configuration from which no toggle(s) can be proposed.\n");
	return WtMCMC_MH_FAILED;
	
      case MH_UNSUCCESSFUL:
	warning("MH MHProposal function failed to find a valid proposal.");
	unsuccessful++;
	if(unsuccessful>taken*MH_QUIT_UNSUCCESSFUL){
	  Rprintf("Too many MH MHProposal function failures.\n");
	  return WtMCMC_MH_FAILED;
	}
      case MH_CONSTRAINT:
	continue;
      }
    }
    
    if(fVerbose>=5){
      Rprintf("MHproposal: ");
      for(unsigned int i=0; i<MHp->ntoggles; i++)
	Rprintf(" (%d, %d)", MHp->toggletail[i], MHp->togglehead[i]);
      Rprintf("\n");
    }

    /* Calculate change statistics,
     remembering that tail -> head */
    WtChangeStats(MHp->ntoggles, MHp->toggletail, MHp->togglehead, MHp->toggleweight, nwp, m);

    /* Always store the proposal for self-tuning. */
    for (unsigned int i = 0; i < m->n_stats; i++){
      prop_networkstatistics[i] += m->workspace[i];
    }

    if(fVerbose>=5){
      Rprintf("Changes: (");
      for(unsigned int i=0; i<m->n_stats; i++)
	Rprintf(" %f ", m->workspace[i]);
      Rprintf(")\n");
    }
    
    /* Calculate the change in the (s-t) %*% W %*% (s-t) due to the proposal. */
    double ip=0;
    for (unsigned int i=0; i<m->n_stats; i++){
     deltainvsig[i]=0.0;
     for (unsigned int j=0; j<m->n_stats; j++){
      deltainvsig[i]+=(m->workspace[j])*invcov[i+(m->n_stats)*j];
     }
     ip+=deltainvsig[i]*((m->workspace[i])+2.0*networkstatistics[i]);
    }
    
    double offsetcontrib = 0;
    for(int i = 0; i < m->n_stats; i++){
        offsetcontrib += (m->workspace[i])*offsets[i];
    }
    
    if(fVerbose>=5){
      Rprintf("log acceptance probability: %f\n", ip - offsetcontrib);
    }
    
    /* if we accept the proposed network */
    if (tau[0]==0? ip - offsetcontrib <= 0 : ip/tau[0] - offsetcontrib <= -log(unif_rand()) ) { 
      if(fVerbose>=5){
	Rprintf("Accepted.\n");
      }

      /* Make proposed toggles (updating timestamps--i.e., for real this time) */
      for(unsigned int i=0; i < MHp->ntoggles; i++){
	Vertex t=MHp->toggletail[i], h=MHp->togglehead[i];
	double w=MHp->toggleweight[i];
	
	if(MHp->discord)
	  for(WtNetwork **nwd=MHp->discord; *nwd!=NULL; nwd++){
	    // This could be speeded up by implementing an "incrementation" function.
	    WtSetEdge(t, h, WtGetEdge(t,  h, *nwd) + w - WtGetEdge(t, h, nwp), *nwd);
	  }

	WtSetEdge(t, h, w, nwp);
      }
      /* record network statistics for posterity */
      Rboolean found = TRUE;
      for (unsigned int i = 0; i < m->n_stats; i++){
	if((networkstatistics[i] += m->workspace[i])!=0) found=FALSE;
      }
      
      taken++;

      if(found)	break;
    }else{
      if(fVerbose>=5){
	Rprintf("Rejected.\n");
      }
    }
  }

  Free(deltainvsig);

  *staken = taken;
  return WtMCMC_OK;
}
