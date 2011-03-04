#include "MCMC.h"

/*****************
 Note on undirected networks:  For j<k, edge {j,k} should be stored
 as (j,k) rather than (k,j).  In other words, only directed networks
 should have (k,j) with k>j.
*****************/

/*****************
 void MCMC_wrapper

 Wrapper for a call from R.
*****************/
void MCMC_wrapper (int *dnumnets, int *nedges,
		   int *heads, int *tails,
                   int *maxpossibleedges,
                   int *dn, int *dflag, int *bipartite, 
                   int *nterms, char **funnames,
                   char **sonames, 
                   char **MHproposaltype, char **MHproposalpackage,
                   double *inputs, double *theta0, int *samplesize, 
                   double *sample, int *burnin, int *interval,  
                   int *newnetworkheads, 
                   int *newnetworktails, 
                   int *fVerbose, 
                   int *attribs, int *maxout, int *maxin, int *minout,
                   int *minin, int *condAllDegExact, int *attriblength, 
                   int *maxedges){
  int directed_flag, hammingterm;
  Vertex n_nodes, nmax, bip, hhead, htail;
  Edge n_networks, nddyads, kedge;
  Network nw[3];
  Model *m;
  ModelTerm *thisterm;
  MHproposal MH;
  
  n_nodes = (Vertex)*dn; 
  n_networks = (Edge)*dnumnets; 
  nmax = (Edge)*maxedges; 
  bip = (Vertex)*bipartite; 
  
  GetRNGstate();  /* R function enabling uniform RNG */
  
  directed_flag = *dflag;

  m=ModelInitialize(*funnames, *sonames, &inputs, *nterms);

  //Rprintf("Init: 0 edges %d yplus edges %d\n",nedges[0],nedges[1]); 
  /* Form the network */
  nw[0]=NetworkInitialize(heads, tails, nedges[0], 
                          n_nodes, directed_flag, bip, 0);
  /* Form the DTERGM network */
  if (nedges[2]>0) {
   heads += nedges[0]+nedges[1];
   tails += nedges[0]+nedges[1];
   nw[2]=NetworkInitialize(heads, tails, nedges[2],
                           n_nodes, directed_flag, bip, 0);
   heads -= nedges[0]+nedges[1];
   tails -= nedges[0]+nedges[1];
  }

  hammingterm=ModelTermHamming (*funnames, *nterms);
  if(hammingterm>0){
/*	     Rprintf("start with setup\n"); */
   Network nwhamming;
   thisterm = m->termarray + hammingterm - 1;
   nddyads = (Edge)(thisterm->inputparams[0]);
   nwhamming=NetworkInitializeD(thisterm->inputparams+1, 
				thisterm->inputparams+1+nddyads, nddyads, 
        n_nodes, directed_flag, bip,0);
   nddyads=0;
   nw[1]=NetworkInitializeD(thisterm->inputparams+1, 
			   thisterm->inputparams+1+nddyads, nddyads,
         n_nodes, directed_flag, bip,0);
/*	     Rprintf("made hw[1]\n"); */
   for (kedge=1; kedge <= nwhamming.nedges; kedge++) {
     FindithEdge(&hhead, &htail, kedge, &nwhamming);
     if(EdgetreeSearch(hhead, htail, nw[0].outedges) == 0){
/*	     Rprintf(" in g0 not g hhead %d htail %d\n",hhead, htail); */
       ToggleEdge(hhead, htail, &nw[1]);
     }
   }
   for (kedge=1; kedge <= nw[0].nedges; kedge++) {
     FindithEdge(&hhead, &htail, kedge, &nw[0]);
     if(EdgetreeSearch(hhead, htail, nwhamming.outedges) == 0){
/*	     Rprintf("not g0  in g hhead %d htail %d\n",hhead, htail); */
       ToggleEdge(hhead, htail, &nw[1]);
     }
   }
/*   Rprintf("Initial number of discordant %d Number of g0 ties %d Number of ties in g %d\n",nw[1].nedges, nwhamming.nedges,nw[0].nedges); */
   NetworkDestroy(&nwhamming);
  }

/* Really this is a formation term */
//	     Rprintf("proposal %s\n",*MHproposaltype); 
//	     Rprintf("proposal %s\n",*MHproposaltype); 
  if(!strncmp(*MHproposaltype,"FormationMLE",12)){
     Rprintf("formation: y0 edges %d yplus edges %d\n",nedges[1],nedges[0]); 
//   Rprintf("proposal %d\n",strncmp(*MHproposaltype,"FormationMLE",12)); 
// nw[1]=NetworkInitialize(heads, tails, nedges[0], 
//                         n_nodes, directed_flag, bip, 0);
  }
  if(!strncmp(*MHproposaltype,"DissolutionMLE",14)){
     Rprintf("dissolution: y0 edges %d yplus edges %d\n",nedges[1],nedges[0]); 
  }
  
  MH_init(&MH,
	  *MHproposaltype, *MHproposalpackage,
	  inputs,
	  *fVerbose,
	  nw, attribs, maxout, maxin, minout, minin,
	  *condAllDegExact, *attriblength);

  MCMCSample (&MH,
	      theta0, sample, *samplesize,
	      *burnin, *interval,
	      hammingterm,
	      *fVerbose, nw, m);

  MH_free(&MH);
        
/* Rprintf("Back! %d %d\n",nw[0].nedges, nmax); */

  /* record new generated network to pass back to R */
  if(nmax>0 && newnetworkheads && newnetworktails)
    newnetworkheads[0]=newnetworktails[0]=EdgeTree2EdgeList(newnetworkheads+1,newnetworktails+1,nw,nmax-1);
  
  ModelDestroy(m);

  NetworkDestroy(nw);
  if (hammingterm > 0)
    NetworkDestroy(&nw[1]);
  if (nedges[2]>0)
    NetworkDestroy(&nw[2]);
  PutRNGstate();  /* Disable RNG before returning */
}


/*********************
 void MCMCSample

 Using the parameters contained in the array theta, obtain the
 network statistics for a sample of size samplesize.  burnin is the
 initial number of Markov chain steps before sampling anything
 and interval is the number of MC steps between successive 
 networks in the sample.  Put all the sampled statistics into
 the networkstatistics array. 
*********************/
void MCMCSample (MHproposal *MHp,
  double *theta, double *networkstatistics, 
  int samplesize, int burnin, 
  int interval, int hammingterm, int fVerbose,
  Network *nwp, Model *m) {
  int staken, tottaken, ptottaken, originterval;
  int i, j, components, diam;
  
  originterval = interval;
  components = diam = 0;
  nwp->duration_info.MCMCtimer=0;
  
/*  if (fVerbose) {
    Rprintf("Simulating %d stats on %ld networks using %s",
             m->n_stats, burnin + samplesize*interval, MHproposaltype);
  } */

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
   MetropolisHastings(MHp, theta, networkstatistics, burnin, &staken,
		      hammingterm, fVerbose, nwp, m);  
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
      MetropolisHastings (MHp, theta, networkstatistics, interval, &staken,
                           hammingterm, fVerbose, nwp, m);
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
      Rprintf("sampler accepted %6.3f%% of %d proposed steps.\n",
      tottaken*100.0/(1.0*originterval*samplesize), originterval*samplesize); 
    }
  }else{
    if (fVerbose){
      Rprintf("sampler accepted %6.3f%% of %d proposed steps.\n",
      staken*100.0/(1.0*burnin), burnin); 
    }
  }
}

/*********************
 void MetropolisHastings

 In this function, theta is a m->n_stats-vector just as in MCMCSample,
 but now networkstatistics is merely another m->n_stats-vector because
 this function merely iterates nsteps times through the Markov
 chain, keeping track of the cumulative change statistics along
 the way, then returns, leaving the updated change statistics in
 the networkstatistics vector.  In other words, this function 
 essentially generates a sample of size one
*********************/
void MetropolisHastings (MHproposal *MHp,
			 double *theta, double *networkstatistics,
			 int nsteps, int *staken,
			 int hammingterm,
       int fVerbose,
			 Network *nwp,
			 Model *m) {
  int step, taken;
  int i;
  double ip, cutoff;
  
  step = taken = 0;
/*  if (fVerbose)
    Rprintf("Now proposing %d MH steps... ", nsteps); */
  while (step < nsteps) {
    MHp->logratio = 0;
    (*(MHp->func))(MHp, nwp); /* Call MH function to propose toggles */
    
    /* Calculate change statistics. */
    ChangeStats(MHp->ntoggles, MHp->togglehead, MHp->toggletail, nwp, m);
      
    /* Calculate inner product */
    for (i=0, ip=0.0; i<m->n_stats; i++){
      ip += theta[i] * m->workspace[i];
    }
    /* The logic is to set exp(cutoff) = exp(ip) * qratio ,
    then let the MH probability equal min{exp(cutoff), 1.0}.
    But we'll do it in log space instead.  */
    cutoff = ip + MHp->logratio;
      
    /* if we accept the proposed network */
    if (cutoff >= 0.0 || log(unif_rand()) < cutoff) { 
      /* Make proposed toggles (updating timestamps--i.e., for real this time) */
      for (i=0; i < MHp->ntoggles; i++){
        ToggleEdgeWithTimestamp(MHp->togglehead[i], MHp->toggletail[i], nwp);
      }
      //    if(!strncmp(MHproposaltype,"FormationMLE",12) |
      //       !strncmp(MHproposaltype,"DissolutionMLE",14) |
      //hammingterm
      if(hammingterm
      ){
        for (i=0; i < MHp->ntoggles; i++){
          Rprintf("Toggle Discord: h %d t %d\n",MHp->togglehead[i],  MHp->toggletail[i]); 
          ToggleEdge(MHp->togglehead[i],  MHp->toggletail[i], &nwp[1]);  /* Toggle the discord for this edge */
        }
      }
      /* record network statistics for posterity */
/*    Rprintf("change stats:");  */
      for (i = 0; i < m->n_stats; i++){
        networkstatistics[i] += m->workspace[i];
/*      Rprintf("%f ", networkstatistics[i]);  */
      }
/*    Rprintf("\n nedges %d\n", nwp->nedges);  */
      taken++;

    }
/*  Catch massive number of edges caused by degeneracy */
/*  if(nwp->nedges > (100000-1000)){step=nsteps;} */
    step++;
    nwp->duration_info.MCMCtimer++;
  }
/*  if (fVerbose)
    Rprintf("%d taken (MCMCtimer=%d)\n", taken, nwp->duration_info.MCMCtimer); */

  *staken = taken;
}

void MCMCPhase12 (int *heads, int *tails, int *dnedges, 
      int *maxpossibleedges,
		  int *dn, int *dflag, int *bipartite, 
		  int *nterms, char **funnames,
		  char **sonames, 
		  char **MHproposaltype, char **MHproposalpackage,
		  double *inputs, 
		  double *theta0, int *samplesize,
		  double *gain, double *meanstats, int *phase1, int *nsub,
		  double *sample, int *burnin, int *interval,  
		  int *newnetworkheads, 
		  int *newnetworktails, 
		  int *fVerbose, 
		  int *attribs, int *maxout, int *maxin, int *minout,
		  int *minin, int *condAllDegExact, int *attriblength, 
		  int *maxedges,
		  int *mheads, int *mtails, int *mdnedges)  {
  int directed_flag, hammingterm, formationterm;
  int nphase1, nsubphases;
  Vertex n_nodes, bip, hhead, htail;
  Edge n_edges, n_medges, nddyads, kedge, nmax;
  Network nw[2];
  Model *m;
  ModelTerm *thisterm;
  MHproposal MH;
  
  nphase1 = *phase1; 
  nsubphases = *nsub;

  n_nodes = (Vertex)*dn; 
  n_edges = (Edge)*dnedges; 
  n_medges = (Edge)*mdnedges; 
  nmax = (Edge)*maxedges; 
  bip = (Vertex)*bipartite; 
  
  GetRNGstate();  /* R function enabling uniform RNG */
  
  directed_flag = *dflag;

  m=ModelInitialize(*funnames, *sonames, &inputs, *nterms);

  /* Form the missing network */
  nw[0]=NetworkInitialize(heads, tails, n_edges,
                          n_nodes, directed_flag, bip,0);

  hammingterm=ModelTermHamming (*funnames, *nterms);
  if(hammingterm>0){
    /*	     Rprintf("start with setup\n"); */
   Network nwhamming;
   thisterm = m->termarray + hammingterm - 1;
   nddyads = (Edge)(thisterm->inputparams[0]);

   nwhamming=NetworkInitializeD(thisterm->inputparams+1,
			       thisterm->inputparams+1+nddyads, nddyads,
             n_nodes, directed_flag, bip, 0);
   nddyads=0;
   nw[1]=NetworkInitializeD(thisterm->inputparams+1,
			   thisterm->inputparams+1+nddyads, nddyads,
         n_nodes, directed_flag, bip,0);
/*	     Rprintf("made hw[1]\n"); */
   for (kedge=1; kedge <= nwhamming.nedges; kedge++) {
     FindithEdge(&hhead, &htail, kedge, &nwhamming);
     if(EdgetreeSearch(hhead, htail, nw[0].outedges) == 0){
/*	     Rprintf(" in g0 not g hhead %d htail %d\n",hhead, htail); */
       ToggleEdge(hhead, htail, &nw[1]);
     }
   }
   for (kedge=1; kedge <= nw[0].nedges; kedge++) {
     FindithEdge(&hhead, &htail, kedge, &nw[0]);
     if(EdgetreeSearch(hhead, htail, nwhamming.outedges) == 0){
/*	     Rprintf("not g0  in g hhead %d htail %d\n",hhead, htail); */
       ToggleEdge(hhead, htail, &nw[1]);
     }
   }
/*   Rprintf("Initial number of discordant %d Number of g0 ties %d Number of ties in g %d\n",nw[1].nedges, nwhamming.nedges,nw[0].nedges); */
   NetworkDestroy(&nwhamming);
  }

/* Really this is a formation term */
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
         n_nodes, directed_flag, bip,0);
/*	     Rprintf("made hw[1]\n"); */
   for (kedge=1; kedge <= nwformation.nedges; kedge++) {
     FindithEdge(&hhead, &htail, kedge, &nwformation);
     if(EdgetreeSearch(hhead, htail, nw[0].outedges) == 0){
/*	     Rprintf(" in g0 not g hhead %d htail %d\n",hhead, htail); */
       ToggleEdge(hhead, htail, &nw[0]);
     }
   }
   for (kedge=1; kedge <= nw[0].nedges; kedge++) {
     FindithEdge(&hhead, &htail, kedge, &nw[0]);
     if(EdgetreeSearch(hhead, htail, nwformation.outedges) == 0){
/*	     Rprintf("not g0  in g hhead %d htail %d\n",hhead, htail); */
       ToggleEdge(hhead, htail, &nw[1]);
     }
   }
/*   Rprintf("Initial number of discordant %d Number of g0 ties %d Number of ties in g %d\n",nw[1].nedges, nwformation.nedges,nw[0].nedges); */
   hammingterm=1;
   NetworkDestroy(&nwformation);
/*   Rprintf("Initial number (discord) from reference %d Number of original %d\n",nw[1].nedges,nw[0].nedges); */
  }
  
  MH_init(&MH,
	  *MHproposaltype, *MHproposalpackage,
	  inputs,
	  *fVerbose,
	  nw, attribs, maxout, maxin, minout, minin,
	  *condAllDegExact, *attriblength);
  
  MCMCSamplePhase12 (&MH,
	      theta0, *gain, meanstats, nphase1, nsubphases, sample, *samplesize,
	      *burnin, *interval,
	      hammingterm,
	      (int)*fVerbose, nw, m);

  MH_free(&MH);
  
  /* record new generated network to pass back to R */
  if(nmax>0 && newnetworkheads && newnetworktails)
    newnetworkheads[0]=newnetworktails[0]=EdgeTree2EdgeList(newnetworkheads+1,newnetworktails+1,nw,nmax);

  ModelDestroy(m);

  NetworkDestroy(nw);
  if (hammingterm > 0  || formationterm > 0)
    NetworkDestroy(&nw[1]);
  PutRNGstate();  /* Disable RNG before returning */
}

/*********************
 void MCMCSamplePhase12

 Using the parameters contained in the array theta, obtain the
 network statistics for a sample of size samplesize.  burnin is the
 initial number of Markov chain steps before sampling anything
 and interval is the number of MC steps between successive 
 networks in the sample.  Put all the sampled statistics into
 the networkstatistics array. 
*********************/
void MCMCSamplePhase12 (MHproposal *MHp,
  double *theta, double gain, double *meanstats, int nphase1, int nsubphases, double *networkstatistics, 
  int samplesize, int burnin, 
  int interval, int hammingterm, int fVerbose,
  Network *nwp, Model *m) {
  int staken, tottaken, ptottaken;
  int i, j, components, diam, iter=0;
  
/*Rprintf("nsubphases %d\n", nsubphases); */

  components = diam = 0;
  nwp->duration_info.MCMCtimer=0;
  
  if (fVerbose)
/*  Rprintf("Total m->n_stats is %i; total samplesize is %d\n", */
    Rprintf("The number of statistics is %i and the total samplesize is %d\n",
             m->n_stats,samplesize);

  /*********************
  networkstatistics are modified in groups of m->n_stats, and they
  reflect the CHANGE in the values of the statistics from the
  original (observed) network.  Thus, when we begin, the initial 
  values of the first group of m->n_stats networkstatistics should 
  all be zero
  *********************/
  double *ubar, *u2bar, *aDdiaginv;
  ubar = (double *)malloc( m->n_stats * sizeof(double));
  u2bar = (double *)malloc( m->n_stats * sizeof(double));
  aDdiaginv = (double *)malloc( m->n_stats * sizeof(double));
  for (j=0; j < m->n_stats; j++){
    networkstatistics[j] = -meanstats[j];
    ubar[j] = 0.0;
    u2bar[j] = 0.0;
  }

  /*********************
   Burn in step.  While we're at it, use burnin statistics to 
   prepare covariance matrix for Mahalanobis distance calculations 
   in subsequent calls to M-H
   *********************/
/*Rprintf("MCMCSampleDyn pre burnin numdissolve %d\n", *numdissolve); */
  
    staken = 0;
    Rprintf("Starting burnin of %d steps\n", burnin);
    MetropolisHastings (MHp, theta,
		  networkstatistics, burnin, &staken,
      hammingterm,
      fVerbose,
		  nwp, m);
    Rprintf("Phase 1: %d steps (interval = %d)\n", nphase1,interval);
    /* Now sample networks */
    for (i=0; i <= nphase1; i++){
      MetropolisHastings (MHp, theta,
		  networkstatistics, interval, &staken,
      hammingterm,
      fVerbose,
		  nwp, m);
      if(i > 0){
       for (j=0; j<m->n_stats; j++){
        ubar[j]  += networkstatistics[j];
        u2bar[j] += networkstatistics[j]*networkstatistics[j];
/*  Rprintf("j %d ubar %f u2bar %f ns %f\n", j,  ubar[j], u2bar[j], */
/*		  networkstatistics[j]); */
       }
      }
    }
    if (fVerbose){
      Rprintf("Returned from Phase 1\n");
      Rprintf("\n gain times inverse variances:\n");
    }
    for (j=0; j<m->n_stats; j++){
      aDdiaginv[j] = u2bar[j]-ubar[j]*ubar[j]/(1.0*nphase1);
      if( aDdiaginv[j] > 0.0){
        aDdiaginv[j] = nphase1*gain/aDdiaginv[j];
      }else{
	aDdiaginv[j]=0.00001;
      }
      if (fVerbose){ Rprintf(" %f", aDdiaginv[j]);}
    }
    if (fVerbose){ Rprintf("\n"); }
  
    staken = 0;
    tottaken = 0;
    ptottaken = 0;
    
    if (fVerbose){
      Rprintf("Phase 2: (samplesize = %d)\n", samplesize);
    }
    /* Now sample networks */
    for (i=1; i < samplesize; i++){
      
      MetropolisHastings (MHp, theta,
		  networkstatistics, interval, &staken,
      hammingterm,
      fVerbose,
		  nwp, m);
    /* Update theta0 */
/*Rprintf("initial:\n"); */
      for (j=0; j<m->n_stats; j++){
        theta[j] -= aDdiaginv[j] * networkstatistics[j];
      }
/*Rprintf("\n"); */
/*    if (fVerbose){ Rprintf("nsubphases %d i %d\n", nsubphases, i); } */
      if (i==(nsubphases)){
	nsubphases = trunc(nsubphases*2.52) + 1;
        if (fVerbose){
	 iter++;
	 Rprintf("End of iteration %d; Updating the number of sub-phases to be %d\n",iter,nsubphases);
	}
        for (j=0; j<m->n_stats; j++){
          aDdiaginv[j] /= 2.0;
          if (fVerbose){Rprintf("theta_%d = %f; change statistic[%d] = %f\n",
		                 j+1, theta[j], j+1, networkstatistics[j]);}
/*        if (fVerbose){ Rprintf(" %f statsmean %f",  theta[j],(networkstatistics[j]-meanstats[j])); } */
        }
        if (fVerbose){ Rprintf("\n"); }
      }
      /* Set current vector of stats equal to previous vector */
      for (j=0; j<m->n_stats; j++){
/*      networkstatistics[j] -= meanstats[j]; */
        networkstatistics[j+m->n_stats] = networkstatistics[j];
      }
      networkstatistics += m->n_stats;
/*      if (fVerbose){ Rprintf("step %d from %d:\n",i, samplesize);} */
      /* This then adds the change statistics to these values */
      tottaken += staken;
#ifdef Win32
      if( ((100*i) % samplesize)==0 && samplesize > 500){
	R_FlushConsole();
    	R_ProcessEvents();
      }
#endif
      if (fVerbose){
        if( ((3*i) % samplesize)==0 && samplesize > 500){
        Rprintf("Sampled %d from Metropolis-Hastings\n", i);}
      }
      
      if( ((3*i) % samplesize)==0 && tottaken == ptottaken){
        ptottaken = tottaken; 
        Rprintf("Warning:  Metropolis-Hastings algorithm has accepted only "
        "%d steps out of a possible %d\n",  ptottaken-tottaken, i); 
      }
/*      Rprintf("Sampled %d from %d\n", i, samplesize); */

    /*********************
    Below is an extremely crude device for letting the user know
    when the chain doesn't accept many of the proposed steps.
    *********************/
/*    if (fVerbose){ */
/*      Rprintf("Metropolis-Hastings accepted %7.3f%% of %d steps.\n", */
/*	      tottaken*100.0/(1.0*interval*samplesize), interval*samplesize);  */
/*    } */
/*  }else{ */
/*    if (fVerbose){ */
/*      Rprintf("Metropolis-Hastings accepted %7.3f%% of %d steps.\n", */
/*	      staken*100.0/(1.0*burnin), burnin);  */
/*    } */
  }
/*  Rprintf("netstats: %d\n", samplesize); */
/*  for (i=0; i < samplesize; i++){ */
/*   for (j=0; j < m->n_stats; j++){ */
/*      Rprintf("%f ", networkstatistics[j+(m->n_stats)*(i)]); */
/*   } */
/*  Rprintf("\n"); */
/*  } */
  if (fVerbose){
    Rprintf("Phase 3: MCMC-Newton-Raphson\n");
  }

  free(ubar);
  free(u2bar);
}

