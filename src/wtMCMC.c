#include "wtMCMC.h"

/*****************
 Note on undirected networks:  For j<k, edge {j,k} should be stored
 as (j,k) rather than (k,j).  In other words, only directed networks
 should have (k,j) with k>j.
*****************/

/*****************
 void WtMCMC_wrapper

 Wrapper for a call from R.

 and don't forget that tail -> head
*****************/
void WtMCMC_wrapper(int *dnumnets, int *nedges,
		    int *tails, int *heads, double *weights, 
		    int *maxpossibleedges,
		    int *dn, int *dflag, int *bipartite, 
		    int *nterms, char **funnames,
		    char **sonames, 
		    char **MHproposaltype, char **MHproposalpackage,
		    double *inputs, double *theta0, int *samplesize, 
		    double *sample, int *burnin, int *interval,  
		    int *newnetworktails, 
		    int *newnetworkheads, 
		    double *newnetworkweights,
		    int *fVerbose, 
		    int *maxedges,
		    int *status){
  int directed_flag;
  Vertex n_nodes, nmax, bip;
  Edge n_networks;
  WtNetwork nw[1];
  WtModel *m;
  WtMHproposal MH;
  
  n_nodes = (Vertex)*dn; 
  n_networks = (Edge)*dnumnets; 
  nmax = (Edge)abs(*maxedges); 
  bip = (Vertex)*bipartite; 
  
  GetRNGstate();  /* R function enabling uniform RNG */
  
  directed_flag = *dflag;

  m=WtModelInitialize(*funnames, *sonames, &inputs, *nterms);

  /* Form the network */
  nw[0]=WtNetworkInitialize(tails, heads, weights, nedges[0], 
			    n_nodes, directed_flag, bip, 0);

  WtMH_init(&MH,
	    *MHproposaltype, *MHproposalpackage,
	    inputs,
	    *fVerbose,
	    nw);

  *status = WtMCMCSample(&MH,
			 theta0, sample, *samplesize,
			 *burnin, *interval,
			 *fVerbose, nmax, nw, m);

  WtMH_free(&MH);
        
/* Rprintf("Back! %d %d\n",nw[0].nedges, nmax); */

  /* record new generated network to pass back to R */
  if(*status == WtMCMC_OK && nmax>0 && newnetworktails && newnetworkheads)
    newnetworktails[0]=newnetworkheads[0]=WtEdgeTree2EdgeList(newnetworktails+1,newnetworkheads+1,newnetworkweights+1,nw,nmax-1);
  
  WtModelDestroy(m);
  WtNetworkDestroy(nw);
  PutRNGstate();  /* Disable RNG before returning */
}


/*********************
 void WtMCMCSample

 Using the parameters contained in the array theta, obtain the
 network statistics for a sample of size samplesize.  burnin is the
 initial number of Markov chain steps before sampling anything
 and interval is the number of MC steps between successive 
 networks in the sample.  Put all the sampled statistics into
 the networkstatistics array. 
*********************/
WtMCMCStatus WtMCMCSample(WtMHproposal *MHp,
			  double *theta, double *networkstatistics, 
			  int samplesize, int burnin, 
			  int interval, int fVerbose, int nmax,
			  WtNetwork *nwp, WtModel *m) {
  int staken, tottaken, ptottaken;
  int i, j;
    
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
   Burn in step.
   *********************/
/*  Catch more edges than we can return */
  if(WtMetropolisHastings(MHp, theta, networkstatistics, burnin, &staken,
			fVerbose, nwp, m)!=WtMCMC_OK)
    return WtMCMC_MH_FAILED;
  if(nmax!=0 && nwp->nedges >= nmax-1){
    return WtMCMC_TOO_MANY_EDGES;
  }
  
/*   if (fVerbose){ 
       Rprintf(".");
     } */
  
  if (samplesize>1){
    staken = 0;
    tottaken = 0;
    ptottaken = 0;
    
    /* Now sample networks */
    for (i=1; i < samplesize; i++){
      /* Set current vector of stats equal to previous vector */
      for (j=0; j<m->n_stats; j++){
        networkstatistics[j+m->n_stats] = networkstatistics[j];
      }
      networkstatistics += m->n_stats;
      /* This then adds the change statistics to these values */
      
      /* Catch massive number of edges caused by degeneracy */
      if(WtMetropolisHastings(MHp, theta, networkstatistics, interval, &staken,
			    fVerbose, nwp, m)!=WtMCMC_OK)
	return WtMCMC_MH_FAILED;
      if(nmax!=0 && nwp->nedges >= nmax-1){
	return WtMCMC_TOO_MANY_EDGES;
      }
      tottaken += staken;

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
      Rprintf("Sampler accepted %6.3f%% of %d proposed steps.\n",
      tottaken*100.0/(1.0*interval*samplesize), interval*samplesize); 
    }
  }else{
    if (fVerbose){
      Rprintf("Sampler accepted %6.3f%% of %d proposed steps.\n",
      staken*100.0/(1.0*burnin), burnin); 
    }
  }
  return WtMCMC_OK;
}

/*********************
 void MetropolisHastings

 In this function, theta is a m->n_stats-vector just as in WtMCMCSample,
 but now networkstatistics is merely another m->n_stats-vector because
 this function merely iterates nsteps times through the Markov
 chain, keeping track of the cumulative change statistics along
 the way, then returns, leaving the updated change statistics in
 the networkstatistics vector.  In other words, this function 
 essentially generates a sample of size one
*********************/
WtMCMCStatus WtMetropolisHastings (WtMHproposal *MHp,
				 double *theta, double *networkstatistics,
				 int nsteps, int *staken,
				 int fVerbose,
				 WtNetwork *nwp,
				 WtModel *m) {
  
  unsigned int taken=0, unsuccessful=0;
/*  if (fVerbose)
    Rprintf("Now proposing %d MH steps... ", nsteps); */
  for(unsigned int step=0; step < nsteps; step++) {
    MHp->logratio = 0;
    (*(MHp->func))(MHp, nwp); /* Call MH function to propose toggles */

    if(MHp->toggletail[0]==MH_FAILED){
      if(MHp->togglehead[0]==MH_UNRECOVERABLE)
	error("Something very bad happened during proposal. Memory has not been deallocated, so restart R soon.");
      if(MHp->togglehead[0]==MH_IMPOSSIBLE){
	Rprintf("MH Proposal function encountered a configuration from which no toggle(s) can be proposed.\n");
	return WtMCMC_MH_FAILED;
      }
      if(MHp->togglehead[0]==MH_UNSUCCESSFUL){
	warning("MH Proposal function failed to find a valid proposal.");
	unsuccessful++;
	if(unsuccessful>taken*MH_QUIT_UNSUCCESSFUL){
	  Rprintf("Too many MH Proposal function failures.\n");
	  return WtMCMC_MH_FAILED;
	}       

	continue;
      }
    }
    
    /* Calculate change statistics,
       remembering that tail -> head */
    WtChangeStats(MHp->ntoggles, MHp->toggletail, MHp->togglehead, MHp->toggleweight, nwp, m);
    
    /* Calculate inner product */
    double ip=0;
    for (unsigned int i=0; i<m->n_stats; i++){
      ip += theta[i] * m->workspace[i];
    }
    /* The logic is to set cutoff = ip+logratio ,
       then let the MH probability equal min{exp(cutoff), 1.0}.
       But we'll do it in log space instead.  */
    double cutoff = ip + MHp->logratio;
    
    /* if we accept the proposed network */
    if (cutoff >= 0.0 || log(unif_rand()) < cutoff) { 
      /* Make proposed toggles (updating timestamps--i.e., for real this time) */
      for(unsigned int i=0; i < MHp->ntoggles; i++){
	WtSetEdge(MHp->toggletail[i], MHp->togglehead[i], MHp->toggleweight[i], nwp);
	
	if(MHp->discord)
	  for(WtNetwork **nwd=MHp->discord; *nwd!=NULL; nwd++){
	    // This could be speeded up by implementing an "incrementation" function.
	    WtSetEdge(MHp->toggletail[i],  MHp->togglehead[i], MHp->toggleweight[i]-WtGetEdge(MHp->toggletail[i],  MHp->togglehead[i], *nwd) + MHp->toggleweight[i]-WtGetEdge(MHp->toggletail[i],  MHp->togglehead[i], nwp), *nwd);
	  }
      }
      /* record network statistics for posterity */
      for (unsigned int i = 0; i < m->n_stats; i++){
	networkstatistics[i] += m->workspace[i];
      }
      taken++;
    }
  }
  
  *staken = taken;
  return WtMCMC_OK;
}

