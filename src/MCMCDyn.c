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
void MCMCDyn_wrapper (double *heads, double *tails, double *dnedges,
                   double *dn, int *dflag, double *bipartite, 
                   int *nterms, char **funnames, char **sonames, 
                   char **MHproposaltype, char **MHproposalpackage,
                   double *inputs, double *theta0, 
                   int *ndynterms, char **dynfunnames, char **dynsonames, 
                   double *dyninputs, 
		   double *samplesize, 
                   double *sample, double *burnin, double *interval,  
                   int *newnetworkhead, int *newnetworktail, 
                   int *numdissolved, int *dissolvedhead, int *dissolvedtail, 
                   int *fVerbose, 
                   double *gamma, int *dyninterval,
                   int *attribs, int *maxout, int *maxin, int *minout,
                   int *minin, int *condAllDegExact, int *attriblength, 
                   double *maxedges,
                   double *mheads, double *mtails, double *mdnedges,
                   int *mdflag)  {
  int i, nextedge, directed_flag, hammingterm, formationterm;
  Vertex v, k, n_nodes, nmax, bip, hhead, htail;
  Edge n_edges, n_medges, nddyads, kedge;
  Network nw[2];
  DegreeBound *bd;
  Model *m, *mdyn;
  ModelTerm *thisterm;
  
  n_nodes = (Vertex)*dn; /* coerce double *dn to type Vertex */
  n_edges = (Vertex)*dnedges; /* coerce double *dnedges to type Vertex */
  n_medges = (Vertex)*mdnedges; /* coerce double *mdnedges to type Vertex */
  nmax = (Vertex)*maxedges; /* coerce double *maxedges to type Vertex */
  bip = (Vertex)*bipartite; /* coerce double *bipartite to type Vertex */
  
  GetRNGstate();  /* R function enabling uniform RNG */
  
  directed_flag = *dflag;

  Edge numdissolve=13;
  Vertex *dissolvehead, *dissolvetail;
  dissolvehead = (Vertex *)malloc(nmax * sizeof(Vertex));
  dissolvetail = (Vertex *)malloc(nmax * sizeof(Vertex));

  for (i = 0; i < nmax; i++){
    newnetworkhead[i] = 0;
    newnetworktail[i] = 0;
    dissolvehead[i] = 0;
    dissolvetail[i] = 0;
  }

  m=ModelInitialize(*funnames, *sonames, inputs, *nterms);

  /* Form the missing network */
  nw[0]=NetworkInitialize(heads, tails, n_edges, n_nodes, directed_flag, bip);
  if (n_medges>0) {
   nw[1]=NetworkInitialize(mheads, mtails, n_medges, n_nodes, directed_flag, bip);
  }

  hammingterm=ModelTermHamming (*funnames, *nterms);
  if(hammingterm>0){
//	     Rprintf("start with setup\n");
   Network nwhamming;
   thisterm = m->termarray + hammingterm - 1;
   nddyads = (Edge)(thisterm->inputparams[0]);
   double *dhead, *dtail;
   dhead = (double *) malloc(sizeof(double) * nddyads);
   dtail = (double *) malloc(sizeof(double) * nddyads);
   for (i=0; i<nddyads; i++){
    dhead[i] = (Vertex)(thisterm->inputparams[1+        i]);
    dtail[i] = (Vertex)(thisterm->inputparams[1+nddyads+i]);
   }
   nwhamming=NetworkInitialize(dhead, dtail, nddyads, n_nodes, directed_flag, bip);
   nddyads=0;
   nw[1]=NetworkInitialize(dhead, dtail, nddyads, n_nodes, directed_flag, bip);
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
   free(dhead);
   free(dtail);
//   Rprintf("Initial number of discordant %d Number of g0 ties %d Number of ties in g %d\n",nw[1].nedges, nwhamming.nedges,nw[0].nedges);
   NetworkDestroy(&nwhamming);
  }

// Really this is a formation term
  formationterm=(*ndynterms);
  if(formationterm>0){
   formationterm=ModelTermDissolve (*dynfunnames, *ndynterms);
   mdyn=ModelInitialize(*dynfunnames, *dynsonames, dyninputs, *ndynterms);
   Network nwformation;
   thisterm = mdyn->termarray + formationterm - 1;
   nddyads = (Edge)(thisterm->inputparams[0]);
   double *dhead, *dtail;
   dhead = (double *) malloc(sizeof(double) * nddyads);
   dtail = (double *) malloc(sizeof(double) * nddyads);
   for (i=0; i<nddyads; i++){
    dhead[i] = (Vertex)(thisterm->inputparams[1+        i]);
    dtail[i] = (Vertex)(thisterm->inputparams[1+nddyads+i]);
   }
   nwformation=NetworkInitialize(dhead, dtail, nddyads, n_nodes, directed_flag, bip);
   nddyads=0;
   nw[1]=NetworkInitialize(dhead, dtail, nddyads, n_nodes, directed_flag, bip);
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
   free(dhead);
   free(dtail);
//   Rprintf("Initial number of discordant %d Number of g0 ties %d Number of ties in g %d\n",nw[1].nedges, nwformation.nedges,nw[0].nedges);
   hammingterm=1;
   NetworkDestroy(&nwformation);
//   Rprintf("Initial number (discord) from reference %d Number of original %d\n",nw[1].nedges,nw[0].nedges);
  }
  
  bd=DegreeBoundInitialize(attribs, maxout, maxin, minout, minin,
			   *condAllDegExact, *attriblength, nw);
  MCMCSampleDyn (*MHproposaltype, *MHproposalpackage,
	      theta0, sample, (long int)*samplesize,
	      (long int)*burnin, (long int)*interval,
	      hammingterm,
	      (int)*fVerbose, gamma, (int)*dyninterval,
	      &numdissolve, dissolvehead, dissolvetail,
	      nw, m, mdyn, bd);
   
// Rprintf("samplesize: %f\n", *samplesize);
// for (i=0; i < (*samplesize)*9; i++){
// 	if(i == 9*trunc(i/9)){Rprintf("\n");}
// Rprintf("%f ", sample[i]);
// }
// Rprintf("\n");

  /* record new generated network to pass back to R */
  nextedge=1;
  if (nw[0].directed_flag) {
   for (v=1; v<=n_nodes; v++) 
    {
      Vertex e;
      for(e = EdgetreeMinimum(nw[0].outedges, v);
	  nw[0].outedges[e].value != 0 && nextedge < nmax;
	  e = EdgetreeSuccessor(nw[0].outedges, e))
	{
          newnetworkhead[nextedge] = v;
          newnetworktail[nextedge] = nw[0].outedges[e].value;
	  nextedge++;
	}
   }
  }else{
   for (v=1; v<=n_nodes; v++) 
    {
      Vertex e;
      for(e = EdgetreeMinimum(nw[0].outedges, v);
	  nw[0].outedges[e].value != 0 && nextedge < nmax;
	  e = EdgetreeSuccessor(nw[0].outedges, e))
	{
          k = nw[0].outedges[e].value;
	  if(v < k){
      newnetworkhead[nextedge] = k;
      newnetworktail[nextedge] = v;
      nextedge++;
	  }else{
      newnetworkhead[nextedge] = v;
      newnetworktail[nextedge] = k;
      nextedge++;
	  }
	}
     }
  }
  newnetworkhead[0]=nextedge;

  *numdissolved=(int)numdissolve;
//Rprintf("numdissolved %d numdissolve %d\n", *numdissolved, numdissolve);
  for (i = 0; i < *numdissolved; i++){
    dissolvedhead[i] = (int)(dissolvehead[i]);
    dissolvedtail[i] = (int)(dissolvetail[i]);
//Rprintf(" %d %d\n", dissolvedhead[i], dissolvedtail[i]);
  }

  ModelDestroy(m);
//  Rprintf("nw edges: %d\n", nw[0].nedges);
  DegreeBoundDestroy(bd);
  NetworkDestroy(nw);
  if (n_medges>0 || hammingterm > 0  || formationterm > 0)
    NetworkDestroy(&nw[1]);
  PutRNGstate();  /* Disable RNG before returning */
}


/*********************
 void MCMCSampleDyn

 Using the parameters contained in the array theta, obtain the
 network statistics for a sample of size samplesize.  burnin is the
 initial number of Markov chain steps before sampling anything
 and interval is the number of MC steps between successive 
 networks in the sample.  Put all the sampled statistics into
 the networkstatistics array. 
*********************/
void MCMCSampleDyn (char *MHproposaltype, char *MHproposalpackage,
  double *theta, double *networkstatistics, 
  long int samplesize, long int burnin, 
  long int interval, int hammingterm, int fVerbose,
  double *gamma, int dyninterval,
  Edge *numdissolve, Vertex *dissolvehead, Vertex *dissolvetail,
  Network *nwp, Model *m, Model *mdyn, DegreeBound *bd) {
  long int staken, tottaken, ptottaken;
  int i, j, components, diam;
  ModelTerm *mtp;
  char *fn, *sn;
  MHproposal MH;
  
  components = diam = 0;
  nwp->duration_info.MCMCtimer=0;
  
  if (fVerbose)
    Rprintf("Total m->n_stats is %i; total samplesize is %d\n",
             m->n_stats,samplesize);

  for (i = 0; MHproposaltype[i] != ' ' && MHproposaltype[i] != 0; i++);
  MHproposaltype[i] = 0;
  /* Extract the required string information from the relevant sources */
  if((fn=(char *)malloc(sizeof(char)*(i+4)))==NULL){
    Rprintf("Error in MCMCSample: Can't allocate %d bytes for fn.\n",
	    sizeof(char)*(i+4));
    exit(0);
  }
  fn[0]='M';
  fn[1]='H';
  fn[2]='_';
  for(j=0;j<i;j++)
    fn[j+3]=MHproposaltype[j];
  fn[i+3]='\0';
  /* fn is now the string 'MH_[name]', where [name] is MHproposaltype */
  for (i = 0; MHproposalpackage[i] != ' ' && MHproposalpackage[i] != 0; i++);
  MHproposalpackage[i] = 0;
  if((sn=(char *)malloc(sizeof(char)*(i+1)))==NULL){
    Rprintf("Error in ModelInitialize: Can't allocate %d bytes for sn.\n",
	    sizeof(char)*(j+1));
    exit(0);
  }
  sn=strncpy(sn,MHproposalpackage,i);
  sn[i]='\0';
  if (fVerbose) 
    Rprintf("MH proposal function is %s from %s package\n",fn,sn);

  /* Search for the MH proposal function pointer */
  MH.func=(void (*)(MHproposal*, DegreeBound*, Network*)) R_FindSymbol(fn,sn,NULL);
  if(MH.func==NULL){
    Rprintf("Error in MCMCSample: could not find function %s in "
	    "namespace for package %s.\n",fn,sn);
    exit(0);
  }      

  /*Clean up by freeing sn and fn*/
  free((void *)fn);
  free((void *)sn);

  MH.ntoggles=0;
  (*(MH.func))(&MH, bd, nwp); /* Call MH proposal function to initialize */
  MH.togglehead = (Vertex *)malloc(MH.ntoggles * sizeof(Vertex));
  MH.toggletail = (Vertex *)malloc(MH.ntoggles * sizeof(Vertex));
  
  /*********************
  networkstatistics are modified in groups of m->n_stats, and they
  reflect the CHANGE in the values of the statistics from the
  original (observed) network.  Thus, when we begin, the initial 
  values of the first group of m->n_stats networkstatistics should 
  all be zero
  *********************/
  for (j=0; j < m->n_stats; j++)
    networkstatistics[j] = 0.0;
  mtp = m->termarray;

  /*********************
   Burn in step.  While we're at it, use burnin statistics to 
   prepare covariance matrix for Mahalanobis distance calculations 
   in subsequent calls to M-H
   *********************/
//Rprintf("MCMCSampleDyn pre burnin numdissolve %d\n", *numdissolve);
  MetropolisHastingsDyn(&MH, theta, networkstatistics, burnin, &staken,
		     hammingterm, fVerbose, gamma, dyninterval,
		     numdissolve, dissolvehead, dissolvetail,
		     nwp, m, mdyn, bd);
//Rprintf("MCMCSampleDyn post burnin numdissolve %d\n", *numdissolve);
  
  if (fVerbose){
    Rprintf("Returned from Metropolis-Hastings burnin\n");
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
      
      MetropolisHastingsDyn (&MH, theta, networkstatistics, interval, &staken,
		  hammingterm, fVerbose, gamma, dyninterval,
		  numdissolve, dissolvehead, dissolvetail,
		  nwp, m, mdyn, bd);
//Rprintf("MCMCSampleDyn loop numdissolve %d\n", *numdissolve);
      tottaken += staken;
      if (fVerbose){
        if( ((3*i) % samplesize)==0 && samplesize > 500){
        Rprintf("Sampled %d from Metropolis-Hastings\n", i);}
      }
      
      if( ((3*i) % samplesize)==0 && tottaken == ptottaken){
        ptottaken = tottaken; 
        Rprintf("Warning:  Metropolis-Hastings algorithm has accepted only "
        "%d steps out of a possible %d\n",  ptottaken-tottaken, i); 
      }
//      Rprintf("Sampled %d from %d\n", i, samplesize);
    }
    /*********************
    Below is an extremely crude device for letting the user know
    when the chain doesn't accept many of the proposed steps.
    *********************/
//    if (fVerbose){
//      Rprintf("Metropolis-Hastings accepted %7.3f%% of %d steps.\n",
//	      tottaken*100.0/(1.0*interval*samplesize), interval*samplesize); 
//    }
//  }else{
//    if (fVerbose){
//      Rprintf("Metropolis-Hastings accepted %7.3f%% of %d steps.\n",
//	      staken*100.0/(1.0*burnin), burnin); 
//    }
  }
//  Rprintf("netstats: %d\n", samplesize);
//  for (i=0; i < samplesize; i++){
//   for (j=0; j < m->n_stats; j++){
//      Rprintf("%f ", networkstatistics[j+(m->n_stats)*(i)]);
//   }
//  Rprintf("\n");
//  }
  free(MH.togglehead);
  free(MH.toggletail);
}

/*********************
 void MetropolisHastingsDyn

 In this function, theta is a m->n_stats-vector just as in MCMCSampleDyn,
 but now networkstatistics is merely another m->n_stats-vector because
 this function merely iterates nsteps times through the Markov
 chain, keeping track of the cumulative change statistics along
 the way, then returns, leaving the updated change statistics in
 the networkstatistics vector.  In other words, this function 
 essentially generates a sample of size one
*********************/
void MetropolisHastingsDyn (MHproposal *MHp,
			 double *theta, double *networkstatistics,
			 long int nsteps, long int *staken,
			 int hammingterm, int fVerbose,
			 double *gamma, int dyninterval, 
			 Edge *numdissolve,
			 Vertex *dissolvehead, Vertex *dissolvetail,
			 Network *nwp,
                         Model *m, Model *mdyn, DegreeBound *bd) {
  Vertex step, dstep;
  int i, curstat=0;
  double *dstats, ip, cutoff;
  ModelTerm *mtp;
//  MHproposal MHdissolve;
//  Vertex dissolventoggles;
  
  Vertex head, tail;
  Edge j, rane, nedges, numdissolved;
  
//  MHdissolve.ntoggles=10000;
//  MHdissolve.togglehead = (Vertex *)malloc(MHdissolve.ntoggles * sizeof(Vertex));
//  MHdissolve.toggletail = (Vertex *)malloc(MHdissolve.ntoggles * sizeof(Vertex));

  step = 0;
  dstep = 0;
  while (dstep < nsteps) {
/*  Dissolve ties */
    nedges=nwp[0].nedges;

//  Rprintf("nterms %d\n", mdyn->n_terms); 
    if (mdyn->n_terms <= 2) { 
     numdissolved = (Edge)(nedges*exp(gamma[0])/(1+exp(gamma[0])));
//    Rprintf("ntoggles %d\n", numdissolved); 
//    Rprintf("gamma[0] %f\n", gamma[0]); 
     if(numdissolved > 0){
     // Sample numdissolved edges without replacement
      for (i = 0; i < nedges; i++){dissolvehead[i] = i;}
      for (j = 0; j < numdissolved; j++) {
	rane = nedges * unif_rand();
	dissolvetail[j] = dissolvehead[rane] + 1;
	dissolvehead[rane] = dissolvehead[--nedges];
      }
      for (i=0; i < numdissolved; i++) {
//       rane = 1 + unif_rand() * nedges;
       FindithEdge(&head, &tail, dissolvetail[i], &nwp[0]);
//       MHdissolve.togglehead[i] = head;
//       MHdissolve.toggletail[i] = tail;
       dissolvehead[i] = head;
       dissolvetail[i] = tail;
      }
     }
    }else{
     numdissolved=0; 
     for (rane=1; rane <= nedges; rane++) {
      FindithEdge(&head, &tail, rane, &nwp[0]);
      dstats = mdyn->workspace;
      mtp = mdyn->termarray;
      for (i=0; i < (mdyn->n_terms-1); i++) {
       /* Calculate change statistics */
       mtp->dstats = dstats;
       (*(mtp->func))(1, &head, &tail, mtp, nwp);
       dstats += (mtp++)->nstats;
      }
// Rprintf("%d: ", rane); 
      for (i=0, ip=0.0; i<(mdyn->n_stats-1); i++){
       ip -= gamma[i] * mdyn->workspace[i];
//  Rprintf("%f ", mdyn->workspace[i]); 
      }
//  Rprintf(" ip=%f\n", ip); 
      if (log(unif_rand()) < ip) { 
       /* Toggles off this edge */
//     MHdissolve.togglehead[MHdissolve.ntoggles] = head;
//     MHdissolve.toggletail[MHdissolve.ntoggles] = tail;
       dissolvehead[numdissolved] = head;
       dissolvetail[numdissolved] = tail;
       numdissolved++; 
      }
     }
     if(numdissolved > 0){numdissolved--;} 
    }

    if(numdissolved > 0){
//  Rprintf("ntoggles=%d\n", numdissolved); 
//      for (i=0; i < numdissolved; i++) {
//Rprintf("%d h %d t %d\n", i, 
//		 dissolvehead[i], dissolvetail[i]); 
//      }
      dstats = m->workspace;
      mtp = m->termarray;
      for (i=0; i < m->n_terms; i++) {
       /* Calculate change statistics */
       mtp->dstats = dstats;
       (*(mtp->func))(numdissolved, dissolvehead, dissolvetail, mtp, nwp);
       dstats += (mtp++)->nstats;
      }
      for (i = 0; i < m->n_stats; i++)
       networkstatistics[i] += m->workspace[i];
    
      for (i=0; i < numdissolved; i++) {
       ToggleEdge(dissolvehead[i], dissolvetail[i], &nwp[0]);
       ToggleEdge(dissolvehead[i], dissolvetail[i], &nwp[1]);
      }
     }

    *numdissolve = numdissolved;
//Rprintf("C numdissolve %d\n", *numdissolve);
//Rprintf("C numdissolved %d\n", numdissolved);

//    Rprintf("dissolve networkstatistics[0] %f\n", networkstatistics[0]); 
//    Rprintf("Made dissolve: start %d end %d\n", nedges, nwp[0].nedges); 

//    MetropolisHastings (&MHdissolve, theta, &networkstatistics, dyninterval, &staken,
//		        hammingterm, fVerbose, nwp, m, bd);

//    for(i=0; i < m->n_stats; i++){
//      m->workspace[i] = 0.0;
//    }

    MHp->ntoggles = 0;
    (*(MHp->func))(MHp, bd, nwp); /* Call MH proposal function to initialize */

    step = 0;
    while (step < dyninterval) {

     MHp->ratio = 1.0;
     (*(MHp->func))(MHp, bd, nwp); /* Call MH function to propose toggles */
     //      Rprintf("Back from proposal; step=%d\n",step);
     mtp = m->termarray;
     dstats = m->workspace;
     curstat = 0;
     
     for (i=0; i < m->n_terms; i++) {
      /* Calculate change statistics */
      mtp->dstats = dstats;
      (*(mtp->func))(MHp->ntoggles, MHp->togglehead, MHp->toggletail, 
                     mtp, nwp);  /* Call d_??? function */
      curstat += (mtp->nstats);
      dstats += (mtp++)->nstats;
     }
      
//  Rprintf("change stats:"); 
     /* Calculate inner product */
     for (i=0, ip=0.0; i<m->n_stats; i++){
      ip += theta[i] * m->workspace[i];
//  Rprintf("%f ", m->workspace[i]); 
     }
//  Rprintf("\n ip %f dedges %f\n", ip, m->workspace[0]); 
     /* The logic is to set exp(cutoff) = exp(ip) * qratio ,
     then let the MH probability equal min{exp(cutoff), 1.0}.
     But we'll do it in log space instead.  */
     cutoff = ip + log(MHp->ratio);
      
     /* if we accept the proposed network */
     if (cutoff >= 0.0 || log(unif_rand()) < cutoff) { 
      /* Make proposed toggles (updating timestamps--i.e., for real this time) */
      for (i=0; i < MHp->ntoggles; i++){
//        ToggleEdgeWithTimestamp(MHp->togglehead[i], MHp->toggletail[i], nwp);
        ToggleEdge(MHp->togglehead[i], MHp->toggletail[i], &nwp[0]);
	ToggleEdge(MHp->togglehead[i],  MHp->toggletail[i], &nwp[1]);  /* Toggle the discord for this edge */
      }
      /* record network statistics for posterity */
      for (i = 0; i < m->n_stats; i++)
        networkstatistics[i] += m->workspace[i];
     }
     step++;
//   nwp->duration_info.MCMCtimer++;
    }

//    Rprintf("End MH: step %d edges %d\n", step, nwp[0].nedges); 
    dstep++;
  }
/*  if (fVerbose)
    Rprintf("%d taken (MCMCtimer=%d)\n", taken, nwp->duration_info.MCMCtimer); */
//    Rprintf("End MH: edges %d discord %d\n", nwp[0].nedges, nwp[1].nedges); 
//    Rprintf("networkstatistics[0] %f\n", networkstatistics[0]); 

//  free(MHdissolve.togglehead);
//  free(MHdissolve.toggletail);

  *staken = dyninterval;
}
