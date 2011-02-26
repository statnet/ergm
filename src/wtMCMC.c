#include "wtMCMC.h"

/*****************
 Note on undirected networks:  For j<k, edge {j,k} should be stored
 as (j,k) rather than (k,j).  In other words, only directed networks
 should have (k,j) with k>j.
*****************/

/*****************
 void WtMCMC_wrapper

 Wrapper for a call from R.
*****************/
void WtMCMC_wrapper (int *dnumnets, int *nedges,
		     int *heads, int *tails, double *weights, 
		     int *maxpossibleedges,
		     int *dn, int *dflag, int *bipartite, 
		     int *nterms, char **funnames,
		     char **sonames, 
		     char **MHproposaltype, char **MHproposalpackage,
		     double *inputs, double *theta0, int *samplesize, 
		     double *sample, int *burnin, int *interval,  
		     int *newnetworkheads, 
		     int *newnetworktails, 
		     double *newnetworkweights,
		     int *fVerbose, 
		     int *maxedges) {
  int directed_flag;
  Vertex n_nodes, nmax, bip;
  Edge n_networks;
  WtNetwork nw[2];
  WtModel *m;
  WtMHproposal MH;
  
  n_nodes = (Vertex)*dn; /* coerce double *dn to type Vertex */
  n_networks = (Edge)*dnumnets; /* coerce double *dnumnets to type Edge */
  nmax = (Edge)*maxedges; /* coerce double *maxedges to type Edge */
  bip = (Vertex)*bipartite; /* coerce double *bipartite to type Vertex */
  
  GetRNGstate();  /* R function enabling uniform RNG */
  
  directed_flag = *dflag;

  m=WtModelInitialize(*funnames, *sonames, &inputs, *nterms);

  /* Form the network */
  nw[0]=WtNetworkInitialize(heads, tails, weights, nedges[0], 
			    n_nodes, directed_flag, bip, 0);
  /* Form the missing network */
  if (nedges[1]>0) {
   heads += nedges[0];
   tails += nedges[0];
   weights += nedges[0];
   nw[1]=WtNetworkInitialize(heads, tails, weights, nedges[1],
                           n_nodes, directed_flag, bip, 0);
   heads -= nedges[0];
   tails -= nedges[0];
   weights -= nedges[0];
  }

  /*  if (fVerbose) {
    Rprintf("Simulating %d stats on %ld networks using %s",
             m->n_stats, burnin + samplesize*interval, WtMHproposaltype);
  } */
  WtMH_init(&MH,
	    *MHproposaltype, *MHproposalpackage,
	    inputs,
	    *fVerbose,
	    nw);

  WtMCMCSample(&MH,
	      theta0, sample, *samplesize,
	      *burnin, *interval,
	      *fVerbose, nw, m);

  WtMH_free(&MH);
/*   int ii;
   double mos=0.0;
   for(ii=0; ii < bd->attrcount; ii++) 
     mos += bd->maxout[ii];
   Rprintf("bd -> attrcount = %d, sum = %f\n", ii, mos); */
        
        
/* Rprintf("Back! %d %d\n",nw[0].nedges, nmax); */

  /* record new generated network to pass back to R */
  if(nmax>0 && newnetworkheads && newnetworktails)
    newnetworkheads[0]=newnetworktails[0]=WtEdgeTree2EdgeList(newnetworkheads+1,newnetworktails+1,newnetworkweights+1,nw,nmax-1);
  
  WtModelDestroy(m);
  WtNetworkDestroy(nw);
  if (nedges[1]>0)
    WtNetworkDestroy(&nw[1]);
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
void WtMCMCSample (WtMHproposal *MHp,
  double *theta, double *networkstatistics, 
  int samplesize, int burnin, 
  int interval, int fVerbose,
  WtNetwork *nwp, WtModel *m) {
  int staken, tottaken, ptottaken, originterval;
  int i, j, components, diam;
  
  originterval = interval;
  components = diam = 0;
  
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
/*  Catch massive number of edges caused by degeneracy */
   if(nwp->nedges > (50000-1000)){burnin=1;}
   WtMetropolisHastings(MHp, theta, networkstatistics, burnin, &staken,
			fVerbose, nwp, m);  
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
      if(nwp->nedges > (50000-1000)){interval=1;}
      WtMetropolisHastings (MHp, theta, networkstatistics, interval, &staken,
                           fVerbose, nwp, m);
      tottaken += staken;

#ifdef Win32
      if( ((100*i) % samplesize)==0 && samplesize > 500){
	R_FlushConsole();
    	R_ProcessEvents();
      }
#endif
      /*if (fVerbose){
        if( ((3*i) % samplesize)==0 && samplesize > 500){
        Rprintf("Sampled %d from Metropolis-Hastings\n", i);}
      }
      if( ((3*i) % samplesize)==0 && tottaken == ptottaken){
        ptottaken = tottaken; 
        Rprintf("Warning:  Metropolis-Hastings algorithm has accepted only "
        "%d steps out of a possible %d\n",  ptottaken-tottaken, i); 
      } 
      if (fVerbose && (i % dotinterval)==0) { 
        Rprintf(".");  
      } */
    }
    /*********************
    Below is an extremely crude device for letting the user know
    when the chain doesn't accept many of the proposed steps.
    *********************/
    if (fVerbose){
      Rprintf("Sampler accepted %6.3f%% of %d proposed steps.\n",
      tottaken*100.0/(1.0*originterval*samplesize), originterval*samplesize); 
    }
  }else{
    if (fVerbose){
      Rprintf("Sampler accepted %6.3f%% of %d proposed steps.\n",
      staken*100.0/(1.0*burnin), burnin); 
    }
  }
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
void WtMetropolisHastings (WtMHproposal *MHp,
			 double *theta, double *networkstatistics,
			 int nsteps, int *staken,
			 int fVerbose,
			 WtNetwork *nwp,
			 WtModel *m) {
  int step, taken;
  int i;
  double ip, cutoff;
  
  step = taken = 0;
  /*  if (fVerbose)
    Rprintf("Now proposing %d MH steps... ", nsteps); */
  while (step < nsteps) {
    MHp->ratio = 1.0;
    (*(MHp->func))(MHp, nwp); /* Call MH function to propose toggles */
    
    // If the proposal failed, skip it.
    if(*MHp->togglehead!=MH_FAILED){
      /* Calculate change statistics. */
      WtChangeStats(MHp->ntoggles, MHp->togglehead, MHp->toggletail, MHp->toggleweight, nwp, m);

      /* Calculate inner product */
      for (i=0, ip=0.0; i<m->n_stats; i++){
	ip += theta[i] * m->workspace[i];
      }
      /* The logic is to set exp(cutoff) = exp(ip) * qratio ,
	 then let the MH probability equal min{exp(cutoff), 1.0}.
	 But we'll do it in log space instead.  */
      cutoff = ip + log(MHp->ratio);
      
      /* if we accept the proposed network */    
      if (cutoff >= 0.0 || log(unif_rand()) < cutoff) { 
	/* Make proposed toggles (updating timestamps--i.e., for real this time) */
	for (i=0; i < MHp->ntoggles; i++){
	  WtSetEdge(MHp->togglehead[i], MHp->toggletail[i], MHp->toggleweight[i], nwp);
	}
	/* record network statistics for posterity */
	for (i = 0; i < m->n_stats; i++){
	  networkstatistics[i] += m->workspace[i];	  
	}
	taken++;
      }
    }else{
      // For the moment, just break.
      if(*MHp->toggletail==MH_IMPOSSIBLE || *MHp->toggletail==MH_UNRECOVERABLE) break;
    }

/*  Catch massive number of edges caused by degeneracy */
/*  if(nwp->nedges > (100000-1000)){step=nsteps;} */
    step++;
  }
/*  if (fVerbose)
    Rprintf("%d taken (MCMCtimer=%d)\n", taken, nwp->duration_info.MCMCtimer); */

  *staken = taken;
}

