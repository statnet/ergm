#include "MCMC.h"
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
void SAN_wrapper (int *heads, int *tails, int *dnedges, 
                   int *maxpossibleedges,
                   int *dn, int *dflag, int *bipartite, 
                   int *nterms, char **funnames,
                   char **sonames, 
                   char **MHproposaltype, char **MHproposalpackage,
                   double *inputs, double *theta0, double *tau, 
		   int *samplesize, 
                   double *sample, int *burnin, int *interval,  
                   int *newnetworkheads, 
                   int *newnetworktails, 
                   double *invcov, 
                   int *fVerbose, 
                   int *attribs, int *maxout, int *maxin, int *minout,
                   int *minin, int *condAllDegExact, int *attriblength, 
                   int *maxedges,
                   int *mheads, int *mtails, int *mdnedges)  {
  int directed_flag, hammingterm, formationterm;
  Vertex n_nodes, nmax, bip, hhead, htail;
  Edge n_edges, n_medges, nddyads, kedge;
  Network nw[2];
  DegreeBound *bd;
  Model *m;
  ModelTerm *thisterm;
  
  n_nodes = (Vertex)*dn; /* coerce double *dn to type Vertex */
  n_edges = (Edge)*dnedges; 
  n_medges = (Edge)*mdnedges;
  nmax = (Edge)*maxedges; 
  bip = (Vertex)*bipartite; /* coerce double *bipartite to type Vertex */
  
  GetRNGstate();  /* R function enabling uniform RNG */
  
  directed_flag = *dflag;


  m=ModelInitialize(*funnames, *sonames, inputs, *nterms);

  /* Form the missing network */
  nw[0]=NetworkInitialize(heads, tails, n_edges,
                          n_nodes, directed_flag, bip, 0);
  if (n_medges>0) {
   nw[1]=NetworkInitialize(mheads, mtails, n_medges,
                          n_nodes, directed_flag, bip, 0);
  }

  hammingterm=ModelTermHamming (*funnames, *nterms);
  if(hammingterm>0){
//	     Rprintf("start with setup\n");
   Network nwhamming;
   thisterm = m->termarray + hammingterm - 1;
   nddyads = (Edge)(thisterm->inputparams[0]);
   nwhamming=NetworkInitializeD(thisterm->inputparams+1,
			       thisterm->inputparams+1+nddyads, nddyads,
             n_nodes, directed_flag, bip, 0);
   nddyads=0;
   nw[1]=NetworkInitializeD(thisterm->inputparams+1,
			   thisterm->inputparams+1+nddyads, nddyads,
         n_nodes, directed_flag, bip, 0);
//	     Rprintf("made hw[1]\n");
   for (kedge=1; kedge <= nwhamming.nedges; kedge++) {
     FindithEdge(&hhead, &htail, kedge, &nwhamming);
     if(EdgetreeSearch(hhead, htail, nw[0].outedges) == 0){
//	     Rprintf(" in g0 not g hhead %d htail %d\n",hhead, htail);
       ToggleEdge(hhead, htail, &nw[1]);
     }
   }
   for (kedge=1; kedge <= nw[0].nedges; kedge++) {
     FindithEdge(&hhead, &htail, kedge, &nw[0]);
     if(EdgetreeSearch(hhead, htail, nwhamming.outedges) == 0){
//	     Rprintf("not g0  in g hhead %d htail %d\n",hhead, htail);
       ToggleEdge(hhead, htail, &nw[1]);
     }
   }
//   Rprintf("Initial number of discordant %d Number of g0 ties %d Number of ties in g %d\n",nw[1].nedges, nwhamming.nedges,nw[0].nedges);
   NetworkDestroy(&nwhamming);
  }

// Really this is a formation term
  formationterm=ModelTermFormation (*funnames, *nterms);
  if(formationterm>0){
   Network nwformation;
   thisterm = m->termarray + formationterm - 1;
   nddyads = (Edge)(thisterm->inputparams[0]);
   nwformation=NetworkInitializeD(thisterm->inputparams+1,
				thisterm->inputparams+1+nddyads, nddyads,
        n_nodes, directed_flag, bip, 0);
   nddyads=0;
   nw[1]=NetworkInitializeD(thisterm->inputparams+1,
			   thisterm->inputparams+1+nddyads, nddyads,
         n_nodes, directed_flag, bip, 0);
//	     Rprintf("made hw[1]\n");
   for (kedge=1; kedge <= nwformation.nedges; kedge++) {
     FindithEdge(&hhead, &htail, kedge, &nwformation);
     if(EdgetreeSearch(hhead, htail, nw[0].outedges) == 0){
//	     Rprintf(" in g0 not g hhead %d htail %d\n",hhead, htail);
       ToggleEdge(hhead, htail, &nw[0]);
     }
   }
   for (kedge=1; kedge <= nw[0].nedges; kedge++) {
     FindithEdge(&hhead, &htail, kedge, &nw[0]);
     if(EdgetreeSearch(hhead, htail, nwformation.outedges) == 0){
//	     Rprintf("not g0  in g hhead %d htail %d\n",hhead, htail);
       ToggleEdge(hhead, htail, &nw[1]);
     }
   }
//   Rprintf("Initial number of discordant %d Number of g0 ties %d Number of ties in g %d\n",nw[1].nedges, nwformation.nedges,nw[0].nedges);
   hammingterm=1;
   NetworkDestroy(&nwformation);
//   Rprintf("Initial number (discord) from reference %d Number of original %d\n",nw[1].nedges,nw[0].nedges);
  }
  
  bd=DegreeBoundInitialize(attribs, maxout, maxin, minout, minin,
			   *condAllDegExact, *attriblength, nw);
  SANSample (*MHproposaltype, *MHproposalpackage,
	      theta0, invcov, tau, sample, (long int)*samplesize,
	      (long int)*burnin, (long int)*interval,
	      hammingterm,
	      (int)*fVerbose, nw, m, bd);
  
  /* record new generated network to pass back to R */
  newnetworkheads[0]=newnetworktails[0]=EdgeTree2EdgeList(newnetworkheads+1,newnetworktails+1,nw,nmax);

  ModelDestroy(m);
  DegreeBoundDestroy(bd);
  NetworkDestroy(nw);
  if (n_medges>0 || hammingterm > 0  || formationterm > 0)
    NetworkDestroy(&nw[1]);
  PutRNGstate();  /* Disable RNG before returning */
}


/*********************
 void SANSample

 Using the parameters contained in the array theta, obtain the
 network statistics for a sample of size samplesize.  burnin is the
 initial number of Markov chain steps before sampling anything
 and interval is the number of MC steps between successive 
 networks in the sample.  Put all the sampled statistics into
 the networkstatistics array. 
*********************/
void SANSample (char *MHproposaltype, char *MHproposalpackage,
  double *theta, double *invcov, double *tau, double *networkstatistics, 
  long int samplesize, long int burnin, 
  long int interval, int hammingterm, int fVerbose,
  Network *nwp, Model *m, DegreeBound *bd) {
  long int staken, tottaken, ptottaken;
  int i, j, components, diam;
  MHproposal MH;
  
  components = diam = 0;
  nwp->duration_info.MCMCtimer=0;
  

  if (fVerbose)
    Rprintf("Total m->n_stats is %i; total samplesize is %d\n",
             m->n_stats,samplesize);

  MH_init(&MH, MHproposaltype, MHproposalpackage, fVerbose, nwp, bd);
  
  /*********************
  networkstatistics are modified in groups of m->n_stats, and they
  reflect the CHANGE in the values of the statistics from the
  original (observed) network.  Thus, when we begin, the initial 
  values of the first group of m->n_stats networkstatistics should 
  all be zero
  *********************/
//  for (j=0; j < m->n_stats; j++)
//    networkstatistics[j] = 0.0;

  /*********************
   Burn in step.  While we're at it, use burnin statistics to 
   prepare covariance matrix for Mahalanobis distance calculations 
   in subsequent calls to M-H
   *********************/
  SANMetropolisHastings(&MH, theta, invcov, tau, networkstatistics, burnin, &staken,
		     hammingterm, fVerbose, nwp, m, bd);
  
  if (fVerbose){
    Rprintf("Returned from SAN Metropolis-Hastings burnin\n");
  }
  
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
      
      SANMetropolisHastings (&MH, theta, invcov, tau, networkstatistics, 
		             interval, &staken,
			     hammingterm, fVerbose, nwp, m, bd);
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
    }
    /*********************
    Below is an extremely crude device for letting the user know
    when the chain doesn't accept many of the proposed steps.
    *********************/
    if (fVerbose){
      Rprintf("SAN Metropolis-Hastings accepted %7.3f%% of %d steps.\n",
	      tottaken*100.0/(1.0*interval*samplesize), interval*samplesize); 
    }
  }else{
    if (fVerbose){
      Rprintf("SAN Metropolis-Hastings accepted %7.3f%% of %d steps.\n",
	      staken*100.0/(1.0*burnin), burnin); 
    }
  }
  MH_free(&MH);
}

/*********************
 void SANMetropolisHastings

 In this function, theta is a m->n_stats-vector just as in SANSample,
 but now networkstatistics is merely another m->n_stats-vector because
 this function merely iterates nsteps times through the Markov
 chain, keeping track of the cumulative change statistics along
 the way, then returns, leaving the updated change statistics in
 the networkstatistics vector.  In other words, this function 
 essentially generates a sample of size one
*********************/
void SANMetropolisHastings (MHproposal *MHp,
			    double *theta, double *invcov, 
			    double *tau, double *networkstatistics,
			    long int nsteps, long int *staken,
			    int hammingterm, int fVerbose,
			    Network *nwp,
			    Model *m, DegreeBound *bd) {
  long int step, taken;
  int i,j;
  double ip,dif;
  double *deltainvsig, *delta;
  deltainvsig = (double *)malloc( m->n_stats * sizeof(double));
  delta = (double *)malloc( m->n_stats * sizeof(double));
  
//  div=0.0;
//    Rprintf("\n");
//for (i=0; i<m->n_stats; i++){
//  div += (networkstatistics[i])*(networkstatistics[i]);
//  Rprintf("i %d %f\n",i,networkstatistics[i]);
//}

  step = taken = 0;
/*  if (fVerbose)
    Rprintf("Now proposing %d MH steps... ", nsteps); */
  while (step < nsteps) {
    MHp->ratio = 1.0;
    (*(MHp->func))(MHp, bd, nwp); /* Call MH function to propose toggles */
    //      Rprintf("Back from proposal; step=%d\n",step);

    ChangeStats(MHp->ntoggles, MHp->togglehead, MHp->toggletail, nwp, m);
      
    dif=0.0;
    ip=0.0;
    /* Calculate inner product */
    for (i=0; i<m->n_stats; i++){
     delta[i]=0.0;
     deltainvsig[i]=0.0;
     for (j=0; j<m->n_stats; j++){
      delta[i]+=networkstatistics[j]*invcov[i+(m->n_stats)*j];
      deltainvsig[i]+=(m->workspace[j])*invcov[i+(m->n_stats)*j];
     }
     ip+=deltainvsig[i]*((m->workspace[i])+2.0*networkstatistics[i]);
     dif+=delta[i]*networkstatistics[i];
    }
//Rprintf("i %d j %d ic %f\n",i,j,invcov[i+(m->n_stats)*j]);
//if(i<=1){Rprintf("i %d %f %f\n",i,div,div*div*dif/(tau[i]*asig2[i]));}
//Rprintf(" ip %f dif %f\n",ip,dif);
//Rprintf("ip %f div %f networkstatistics[0] %f networkstatistics[1] %f\n",
//	 ip,div,networkstatistics[0],networkstatistics[1]);
//  Rprintf("ip %f m->workspace[i] %f ns %f asig2[0] %f div %f \n",ip, m->workspace[0],networkstatistics[0],asig2[0],div);
// Rprintf("step %d tau[0] %f tau[1] %f div %f \n",step, tau[0],tau[1],div);
      
    /* if we accept the proposed network */
    if (ip <= 0.0) { 
//  if (ip <= 0.0 || (ip/dif) < 0.001) { 
//  if (div > 0.0 && (ip < 0.0 || unif_rand() < 0.01)) { 
// if (ip <= 0.0 || (ip/dif) < (nsteps-step)*0.001*tau[0]/(1.0*nsteps)) { 
//  if (ip > exp(theta[0])*(m->n_stats)*unif_rand()/(1.0+exp(theta[0])) { 
//  if (ip > tau[0]*(m->n_stats)*unif_rand()) { 
      /* Make proposed toggles (updating timestamps--i.e., for real this time) */
      for (i=0; i < MHp->ntoggles; i++){
        ToggleEdgeWithTimestamp(MHp->togglehead[i], MHp->toggletail[i], nwp);
        if(hammingterm){
	 ToggleEdge(MHp->togglehead[i],  MHp->toggletail[i], &nwp[1]);  /* Toggle the discord for this edge */
	}
      }
      /* record network statistics for posterity */
      for (i = 0; i < m->n_stats; i++){
        networkstatistics[i] += m->workspace[i];
      }
//  div=0.0;
//    for (i=0; i<m->n_stats; i++){
//      div += (networkstatistics[i])*(networkstatistics[i]);
//    }
      taken++;
    }
    step++;
    nwp->duration_info.MCMCtimer++;
  }
/*  if (fVerbose)
    Rprintf("%d taken (MCMCtimer=%d)\n", taken, nwp->duration_info.MCMCtimer); */

  free(deltainvsig);
  free(delta);

  *staken = taken;
}
