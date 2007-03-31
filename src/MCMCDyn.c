#include "MCMC.h"
#include "MCMCDyn.h"

/*****************
 Note on undirected networks:  For j<k, edge {j,k} should be stored
 as (j,k) rather than (k,j).  In other words, only directed networks
 should have (k,j) with k>j.
*****************/

/*****************
 void MCMC_wrapper

 Wrapper for a call from R.
*****************/
void MCMCDyn_wrapper (int *order_code, double *heads, double *tails, double *dnedges,
                   double *dn, int *dflag, double *bipartite, 
                   int *nterms, char **funnames, char **sonames, 
                   char **MHproposaltype, char **MHproposalpackage,
                   double *inputs, double *theta0, 
                   int *ndynterms, char **dynfunnames, char **dynsonames, 
                   double *dyninputs, 
		   double *samplesize, 
                   double *sample, double *burnin, double *interval,  
                   int *newnetworkhead, int *newnetworktail, 
                   int *diffnetworktime, int *diffnetworkhead, int *diffnetworktail, 
                   int *dissnetworktime, int *dissnetworkhead, int *dissnetworktail, 
                   int *fVerbose, 
                   double *gamma, int *dyninterval,
                   int *attribs, int *maxout, int *maxin, int *minout,
                   int *minin, int *condAllDegExact, int *attriblength, 
                   double *maxedges,
                   double *mheads, double *mtails, double *mdnedges,
                   int *mdflag)  {
  int i, nextedge, directed_flag, hammingterm, formationterm;
  Vertex v, k, n_nodes, bip, hhead, htail;
  Edge n_edges, n_medges, nddyads, kedge, nmax;
  Network nw[2];
  DegreeBound *bd;
  Model *m, *mdyn;
  ModelTerm *thisterm;
  DynamOrder order;

  switch(*order_code){
  case 1: order=DissThenForm; break;
  case 2: order=DissAndForm; break;
  default:
    error("Unsupported dynamic model code %d.", order_code);
  }
  
  n_nodes = (Vertex)*dn; /* coerce double *dn to type Vertex */
  n_edges = (Vertex)*dnedges; /* coerce double *dnedges to type Vertex */
  n_medges = (Vertex)*mdnedges; /* coerce double *mdnedges to type Vertex */
  nmax = (Edge)*maxedges; /* coerce double *maxedges to type Edge */
  bip = (Vertex)*bipartite; /* coerce double *bipartite to type Vertex */
  
  GetRNGstate();  /* R function enabling uniform RNG */
  
  directed_flag = *dflag;

  Vertex *dissolvetime, *dissolvehead, *dissolvetail;
  dissolvetime = (Vertex *)malloc(nmax * sizeof(Vertex));
  dissolvehead = (Vertex *)malloc(nmax * sizeof(Vertex));
  dissolvetail = (Vertex *)malloc(nmax * sizeof(Vertex));

  Vertex *disstime, *disshead, *disstail;
  disstime = (Vertex *)malloc(nmax * sizeof(Vertex));
  disshead = (Vertex *)malloc(nmax * sizeof(Vertex));
  disstail = (Vertex *)malloc(nmax * sizeof(Vertex));

  Vertex *difftime, *diffhead, *difftail;
  difftime = (Vertex *)malloc(nmax * sizeof(Vertex));
  diffhead = (Vertex *)malloc(nmax * sizeof(Vertex));
  difftail = (Vertex *)malloc(nmax * sizeof(Vertex));

  for (i = 0; i < nmax; i++){
    newnetworkhead[i] = 0;
    newnetworktail[i] = 0;
    diffnetworktime[i] = 0;
    diffnetworkhead[i] = 0;
    diffnetworktail[i] = 0;
    dissolvetime[i] = 0;
    dissolvehead[i] = 0;
    dissolvetail[i] = 0;
    difftime[i] = 0;
    diffhead[i] = 0;
    difftail[i] = 0;
    disstime[i] = 0;
    disshead[i] = 0;
    disstail[i] = 0;
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

// Really this is a dissolve/formation term
  formationterm=(*ndynterms);
  if(formationterm>0){
   formationterm=ModelTermDissolve (*dynfunnames, *ndynterms);
//   Rprintf("*ndynterms %d formationterm %d\n",*ndynterms, formationterm);
   mdyn=ModelInitialize(*dynfunnames, *dynsonames, dyninputs, *ndynterms);
   Network nwformation;
   thisterm = mdyn->termarray + formationterm - 1;
   nddyads = (Edge)(thisterm->inputparams[0]);
//   Rprintf("nddyads %d\n",nddyads);
   double *dhead, *dtail;
   dhead = (double *) malloc(sizeof(double) * nddyads);
   dtail = (double *) malloc(sizeof(double) * nddyads);
   for (i=0; i<nddyads; i++){
    dhead[i] = (Vertex)(thisterm->attrib[        i]);
    dtail[i] = (Vertex)(thisterm->attrib[nddyads+i]);
   }
//    dhead[i] = (Vertex)(thisterm->inputparams[1+        i]);
//    dtail[i] = (Vertex)(thisterm->inputparams[1+nddyads+i]);
   nwformation=NetworkInitialize(dhead, dtail, nddyads, n_nodes, directed_flag, bip);
   nddyads=0;
   nw[1]=NetworkInitialize(dhead, dtail, nddyads, n_nodes, directed_flag, bip);
//	     Rprintf("made hw[1]\n");
   for (kedge=1; kedge <= nwformation.nedges; kedge++) {
     FindithEdge(&hhead, &htail, kedge, &nwformation);
     if(EdgetreeSearch(hhead, htail, nw[0].outedges) == 0){
//	     Rprintf(" in g0 not g hhead %d htail %d\n",hhead, htail);
       ToggleEdge(hhead, htail, &nw[1]);
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
  MCMCSampleDyn (order, *MHproposaltype, *MHproposalpackage,
	      theta0, sample, (long int)*samplesize,
	      (long int)*burnin, (long int)*interval,
	      hammingterm,
	      (int)*fVerbose, gamma, (int)*dyninterval,
	      &nmax,
	      dissolvetime, dissolvehead, dissolvetail,
	      disstime, disshead, disstail,
	      difftime, diffhead, difftail,
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

  diffnetworktime[0] = (int)(diffhead[0]);
  diffnetworkhead[0] = (int)(diffhead[0]);
  diffnetworktail[0] = (int)(diffhead[0]);
//Rprintf("numdissolved %d numdissolve %d\n", *numdissolved, numdissolve);
  for (i = 1; i <= diffnetworkhead[0]; i++){
    diffnetworktime[i] = (int)(difftime[i]);
    diffnetworkhead[i] = (int)(diffhead[i]);
    diffnetworktail[i] = (int)(difftail[i]);
//Rprintf(" %d %d\n", dissolvedhead[i], dissolvedtail[i]);
  }

  dissnetworktime[0] = (int)(disshead[0]);
  dissnetworkhead[0] = (int)(disshead[0]);
  dissnetworktail[0] = (int)(disshead[0]);
//Rprintf("numdissolved %d numdissolve %d\n", *numdissolved, numdissolve);
  for (i = 1; i <= dissnetworkhead[0]; i++){
    dissnetworktime[i] = (int)(disstime[i]);
    dissnetworkhead[i] = (int)(disshead[i]);
    dissnetworktail[i] = (int)(disstail[i]);
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
void MCMCSampleDyn (DynamOrder order, char *MHproposaltype, char *MHproposalpackage,
  double *theta, double *networkstatistics, 
  long int samplesize, long int burnin, 
  long int interval, int hammingterm, int fVerbose,
  double *gamma, int dyninterval,
  Edge *nmax,
  Vertex *dissolvetime, Vertex *dissolvehead, Vertex *dissolvetail,
  Vertex *disstime, Vertex *disshead, Vertex *disstail,
  Vertex *difftime, Vertex *diffhead, Vertex *difftail,
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
    error("Error in MCMCSample: Can't allocate %d bytes for fn. Memory has not been deallocated, so restart R sometime soon.\n",
	    sizeof(char)*(i+4));
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
    error("Error in ModelInitialize: Can't allocate %d bytes for sn. Memory has not been deallocated, so restart R sometime soon.\n",
	  sizeof(char)*(j+1));
  }
  sn=strncpy(sn,MHproposalpackage,i);
  sn[i]='\0';
  if (fVerbose) 
    Rprintf("MH proposal function is %s from %s package\n",fn,sn);

  /* Search for the MH proposal function pointer */
  MH.func=(void (*)(MHproposal*, DegreeBound*, Network*)) R_FindSymbol(fn,sn,NULL);
  if(MH.func==NULL){
    error("Error in MCMCSample: could not find function %s in "
	  "namespace for package %s."
	  "Memory has not been deallocated, so restart R sometime soon.\n",fn,sn);
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
  MetropolisHastingsDyn(order, &MH, theta, networkstatistics, burnin, &staken,
		     hammingterm, fVerbose, gamma, dyninterval,
		     nmax,
		     dissolvetime, dissolvehead, dissolvetail,
		     disstime, disshead, disstail,
		     difftime, diffhead, difftail,
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
      
      MetropolisHastingsDyn (order, &MH, theta, networkstatistics, interval, &staken,
		  hammingterm, fVerbose, gamma, dyninterval,
		  nmax,
		  dissolvetime, dissolvehead, dissolvetail,
		  disstime, disshead, disstail,
		  difftime, diffhead, difftail,
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

/* Helper function to delete all the edges of a network */
R_INLINE void MetropolisHastingsDyn_null_network(Network *nwp){
  Vertex head, tail;
  while (nwp->nedges>0) { // While a network has edges,
    FindithEdge(&head, &tail, 1, nwp); // find one,
    ToggleEdge(head, tail, nwp); // and delete it.
  }
}

/* Helper function to select ties to be dissolved */
R_INLINE int MetropolisHastingsDyn_choose_dissolved(Edge *nmax, double *gamma,
						    Vertex *dissolvehead, Vertex *dissolvetail,
						    Network *nwp, Model *mdyn){
  unsigned int i;
  Vertex head, tail;
  double *dstats, ip;
  Edge j, rane, nedges, numdissolved;
  ModelTerm *mtp;
  nedges=nwp[0].nedges;
  
  //  Rprintf("nterms %d\n", mdyn->n_terms); 
//  Rprintf("choose *nmax %d\n", *nmax); 
  if (mdyn->n_terms <= 2) { 
    //   fast code for edges-only dissolve model
    //   The next line for a fixed number dissolved
    //   numdissolved = (Edge)(nedges*exp(gamma[0])/(1+exp(gamma[0])));
    //   The next line for a random (Binomial) number dissolved
    while(
      (numdissolved = rbinom((double) nedges, exp(gamma[0])/(1+exp(gamma[0]))))>*nmax);
    //    Rprintf("numdissolved %d\n", numdissolved); 
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
	rane = (int)(dissolvetail[i]);
	FindithEdge(&head, &tail, rane, nwp);
	dissolvehead[i] = head;
	dissolvetail[i] = tail;
      }
    }
  }else{
    // Slow code for more complex dissolve models 
    numdissolved=0; 
    for (rane=1; rane <= nedges; rane++) {
      FindithEdge(&head, &tail, rane, nwp);
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
      if (numdissolved >= *nmax) { Rprintf(" numdissolved=%d\n", numdissolved);} 
      if (log(unif_rand()) < ip && numdissolved < *nmax) { 
	dissolvehead[numdissolved] = head;
	dissolvetail[numdissolved] = tail;
	numdissolved++; 
      }
    }
    if(numdissolved > 0){numdissolved--;} 
  }
  return(numdissolved);
}

/* Helper function to commit the ties to be dissolved */
R_INLINE void MetropolisHastingsDyn_commit_dissolve(double *networkstatistics,
						    double *gamma,
						    Vertex *dissolvehead, Vertex *dissolvetail,
						    Vertex *disstime, Vertex *disshead, Vertex *disstail,
						    Network *nwp,
						    Model *m, Model *mdyn,
						    Vertex dstep, Edge *nmax, Edge numdissolved, Edge *nextdissedge){
  double *dstats;
  ModelTerm *mtp;
  unsigned int i;
//  Rprintf("commit *nmax %d\n", *nmax); 
  //  Rprintf("numdissolved=%d\n", numdissolved); 
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
    for (i=0; i < numdissolved && *nextdissedge < *nmax; i++){
    ToggleEdge(dissolvehead[i], dissolvetail[i], &nwp[0]);
    ToggleEdge(dissolvehead[i], dissolvetail[i], &nwp[1]);
    disstime[*nextdissedge] = dstep;
    disshead[*nextdissedge] = dissolvehead[i];
    disstail[*nextdissedge] = dissolvetail[i];
    (*nextdissedge)++;
  }
}

/* Helper function to run the formation side of the process */
R_INLINE void MetropolisHastingsDyn_form(MHproposal *MHp,
					 double *theta, double *networkstatistics,
					 int dyninterval, 
					 Network *nwp,
					 Model *m, DegreeBound *bd){
  Vertex step;
  double *dstats, cutoff, ip;
  unsigned int i;
  int curstat;
  ModelTerm *mtp;
  
  MHp->ntoggles = 0;
  (*(MHp->func))(MHp, bd, nwp); /* Call MH proposal function to initialize */
  
  for(step = 0; step < dyninterval; step++) {
    
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
	  ToggleEdge(MHp->togglehead[i], MHp->toggletail[i], &nwp[1]);  /* Toggle the discord for this edge */
	}
	/* record network statistics for posterity */
	for (i = 0; i < m->n_stats; i++)
	  networkstatistics[i] += m->workspace[i];
    }
    //   nwp->duration_info.MCMCtimer++;
  }
}

/* Helper function to record new generated network DIFFERENCES to pass back to R */
R_INLINE void MetropolisHastingsDyn_record_diff(Edge *nmax,
						Vertex *difftime, Vertex *diffhead, Vertex *difftail,
						Network *nwp, 
						Vertex dstep, Edge *nextdiffedge){
  Vertex head, v;

  if(nwp->directed_flag) {
    for(v=1; v<=nwp->nnodes; v++){
      Vertex e;
      for(e = EdgetreeMinimum(nwp->outedges, v);
	  nwp->outedges[e].value != 0 && *nextdiffedge < *nmax;
	  e = EdgetreeSuccessor(nwp->outedges, e)){
	difftime[*nextdiffedge] = dstep;
	diffhead[*nextdiffedge] = v;
	difftail[*nextdiffedge] = nwp->outedges[e].value;
	(*nextdiffedge)++;
      }
    }
  }else{
    for(v=1; v<=nwp->nnodes; v++){
      Vertex e;
      for(e = EdgetreeMinimum(nwp->outedges, v);
	  nwp->outedges[e].value != 0 && *nextdiffedge < *nmax;
	  e = EdgetreeSuccessor(nwp->outedges, e)){
	head = nwp->outedges[e].value;
	difftime[*nextdiffedge] = dstep;
	if(v < head){
	  diffhead[*nextdiffedge] = head;
	  difftail[*nextdiffedge] = v;
	  (*nextdiffedge)++;
	}else{
	  diffhead[*nextdiffedge] = v;
	  difftail[*nextdiffedge] = head;
	  (*nextdiffedge)++;
	}
      }
    }
  }
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
 // NOTE: To fit a simultaneous model, move the call to
 // MetropolisHastingsDyn_form to between the calls to
 // _choose_dissolved and _commit_dissolve. -- PK
*********************/
void MetropolisHastingsDyn(DynamOrder order, MHproposal *MHp,
			   double *theta, double *networkstatistics,
			   long int nsteps, long int *staken,
			   int hammingterm, int fVerbose,
			   double *gamma, int dyninterval, 
			   Edge *nmax,
			   Vertex *dissolvetime, Vertex *dissolvehead, Vertex *dissolvetail,
			   Vertex *disstime, Vertex *disshead, Vertex *disstail,
			   Vertex *difftime, Vertex *diffhead, Vertex *difftail,
			   Network *nwp,
			   Model *m, Model *mdyn, DegreeBound *bd) {
  unsigned int dstep;
  Edge nextdiffedge, nextdissedge;

  Edge numdissolved;
  
  nextdiffedge=1;
  nextdissedge=1;
  
//  Rprintf("Dyn *nmax %d nmax %d\n", *nmax, nmax); 

  for(dstep = 0; dstep < nsteps; dstep++){
    switch(order){
    case DissThenForm:
      /* Reset reference discord to null when starting to dissolve*/
      if(dstep > 0){
	MetropolisHastingsDyn_null_network(&nwp[1]);
      }
      /* Choose ties to dissolve. */
      numdissolved = MetropolisHastingsDyn_choose_dissolved(nmax, gamma,
							    dissolvehead, dissolvetail,
							    nwp, mdyn);
      /* If any were, commit the dissolution. */
      if(numdissolved > 0){
	MetropolisHastingsDyn_commit_dissolve(networkstatistics,
					      gamma,
					      dissolvehead, dissolvetail,
					      disstime, disshead, disstail,
					      nwp, m, mdyn,
					      dstep, nmax, numdissolved, &nextdissedge);
      }
      /* Run the formation process. */
      MetropolisHastingsDyn_form(MHp,
				 theta, networkstatistics,
				 dyninterval,
				 nwp,
				 m, bd);
      break;
    case DissAndForm:
       /* Reset reference discord to null when starting to dissolve*/
      if(dstep > 0){
	MetropolisHastingsDyn_null_network(&nwp[1]);
      }
      /* Choose ties to dissolve. */
      numdissolved = MetropolisHastingsDyn_choose_dissolved(nmax, gamma,
							    dissolvehead, dissolvetail,
							    nwp, mdyn);
      /* Run the formation process. */
      MetropolisHastingsDyn_form(MHp,
				 theta, networkstatistics,
				 dyninterval,
				 nwp,
				 m, bd);
      /* Commit the dissolution. */
      if(numdissolved > 0){
	MetropolisHastingsDyn_commit_dissolve(networkstatistics,
					      gamma,
					      dissolvehead, dissolvetail,
					      disstime, disshead, disstail,
					      nwp, m, mdyn,
					      dstep, nmax, numdissolved, &nextdissedge);
      }
      break;
    default: 
      error("Unsupported dynamic model code %d. Memory has not been deallocated so restart R sometime soon.",order);
    }
    /* Record toggled dyads. */
    MetropolisHastingsDyn_record_diff(nmax,
				      difftime, diffhead, difftail,
				      &nwp[1],
				      dstep, &nextdiffedge);
  }

//Rprintf("record diff *nmax %d *nextdiffedge %d\n", *nmax, nextdiffedge); 
  difftime[0]=nextdiffedge;
  diffhead[0]=nextdiffedge;
  difftail[0]=nextdiffedge;
  disstime[0]=nextdissedge;
  disshead[0]=nextdissedge;
  disstail[0]=nextdissedge;
  /*  if (fVerbose)
      Rprintf("%d taken (MCMCtimer=%d)\n", taken, nwp->duration_info.MCMCtimer); */
  //    Rprintf("End MH: edges %d discord %d\n", nwp[0].nedges, nwp[1].nedges); 
  //    Rprintf("networkstatistics[0] %f\n", networkstatistics[0]); 
  
  *staken = dyninterval;
}

void MCMCDynPhase2 (int *order_code, double *heads, double *tails, double *dnedges,
                   double *dn, int *dflag, double *bipartite, 
                   int *nterms, char **funnames, char **sonames, 
                   char **MHproposaltype, char **MHproposalpackage,
                   double *inputs, 
		   double *theta0, double *gain, double *meanstats, 
		   int *phase1, int *nsub,
                   int *ndynterms, char **dynfunnames, char **dynsonames, 
                   double *dyninputs, 
		   double *samplesize, 
                   double *sample, double *burnin, double *interval,  
                   int *newnetworkhead, int *newnetworktail, 
                   int *diffnetworktime, int *diffnetworkhead, int *diffnetworktail, 
                   int *dissnetworktime, int *dissnetworkhead, int *dissnetworktail, 
                   int *fVerbose, 
                   double *gamma, int *dyninterval,
                   int *attribs, int *maxout, int *maxin, int *minout,
                   int *minin, int *condAllDegExact, int *attriblength, 
                   double *maxedges,
                   double *mheads, double *mtails, double *mdnedges,
                   int *mdflag)  {
  int i, nextedge, directed_flag, hammingterm, formationterm;
  int nphase1, nsubphases;
  Vertex v, k, n_nodes, bip, hhead, htail;
  Edge n_edges, n_medges, nddyads, kedge, nmax;
  Network nw[2];
  DegreeBound *bd;
  Model *m, *mdyn;
  ModelTerm *thisterm;
  DynamOrder order;

  switch(*order_code){
  case 1: order=DissThenForm; break;
  case 2: order=DissAndForm; break;
  default:
    error("Unsupported dynamic model code %d.", order_code);
  }
  
  nphase1 = (int)*phase1; /* coerce double *dn to type Vertex */
  nsubphases = (int)*nsub; /* coerce double *dn to type Vertex */
  n_nodes = (Vertex)*dn; /* coerce double *dn to type Vertex */
  n_edges = (Vertex)*dnedges; /* coerce double *dnedges to type Vertex */
  n_medges = (Vertex)*mdnedges; /* coerce double *mdnedges to type Vertex */
  nmax = (Edge)*maxedges; /* coerce double *maxedges to type Edge */
  bip = (Vertex)*bipartite; /* coerce double *bipartite to type Vertex */
  
  GetRNGstate();  /* R function enabling uniform RNG */
  
  directed_flag = *dflag;

  Vertex *dissolvetime, *dissolvehead, *dissolvetail;
  dissolvetime = (Vertex *)malloc(nmax * sizeof(Vertex));
  dissolvehead = (Vertex *)malloc(nmax * sizeof(Vertex));
  dissolvetail = (Vertex *)malloc(nmax * sizeof(Vertex));

  Vertex *disstime, *disshead, *disstail;
  disstime = (Vertex *)malloc(nmax * sizeof(Vertex));
  disshead = (Vertex *)malloc(nmax * sizeof(Vertex));
  disstail = (Vertex *)malloc(nmax * sizeof(Vertex));

  Vertex *difftime, *diffhead, *difftail;
  difftime = (Vertex *)malloc(nmax * sizeof(Vertex));
  diffhead = (Vertex *)malloc(nmax * sizeof(Vertex));
  difftail = (Vertex *)malloc(nmax * sizeof(Vertex));

  for (i = 0; i < nmax; i++){
    newnetworkhead[i] = 0;
    newnetworktail[i] = 0;
    diffnetworktime[i] = 0;
    diffnetworkhead[i] = 0;
    diffnetworktail[i] = 0;
    dissolvetime[i] = 0;
    dissolvehead[i] = 0;
    dissolvetail[i] = 0;
    difftime[i] = 0;
    diffhead[i] = 0;
    difftail[i] = 0;
    disstime[i] = 0;
    disshead[i] = 0;
    disstail[i] = 0;
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

// Really this is a dissolve/formation term
  formationterm=(*ndynterms);
  if(formationterm>0){
   formationterm=ModelTermDissolve (*dynfunnames, *ndynterms);
//   Rprintf("*ndynterms %d formationterm %d\n",*ndynterms, formationterm);
   mdyn=ModelInitialize(*dynfunnames, *dynsonames, dyninputs, *ndynterms);
   Network nwformation;
   thisterm = mdyn->termarray + formationterm - 1;
   nddyads = (Edge)(thisterm->inputparams[0]);
//   Rprintf("nddyads %d\n",nddyads);
   double *dhead, *dtail;
   dhead = (double *) malloc(sizeof(double) * nddyads);
   dtail = (double *) malloc(sizeof(double) * nddyads);
   for (i=0; i<nddyads; i++){
    dhead[i] = (Vertex)(thisterm->attrib[        i]);
    dtail[i] = (Vertex)(thisterm->attrib[nddyads+i]);
   }
//    dhead[i] = (Vertex)(thisterm->inputparams[1+        i]);
//    dtail[i] = (Vertex)(thisterm->inputparams[1+nddyads+i]);
   nwformation=NetworkInitialize(dhead, dtail, nddyads, n_nodes, directed_flag, bip);
   nddyads=0;
   nw[1]=NetworkInitialize(dhead, dtail, nddyads, n_nodes, directed_flag, bip);
//	     Rprintf("made hw[1]\n");
   for (kedge=1; kedge <= nwformation.nedges; kedge++) {
     FindithEdge(&hhead, &htail, kedge, &nwformation);
     if(EdgetreeSearch(hhead, htail, nw[0].outedges) == 0){
//	     Rprintf(" in g0 not g hhead %d htail %d\n",hhead, htail);
       ToggleEdge(hhead, htail, &nw[1]);
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
  
//Rprintf("nsubphases %d\n", nsubphases);

  bd=DegreeBoundInitialize(attribs, maxout, maxin, minout, minin,
			   *condAllDegExact, *attriblength, nw);
  MCMCSampleDynPhase2 (order, *MHproposaltype, *MHproposalpackage,
	      theta0, *gain, meanstats, nphase1, nsubphases, sample, (long int)*samplesize,
	      (long int)*burnin, (long int)*interval,
	      hammingterm,
	      (int)*fVerbose, gamma, (int)*dyninterval,
	      &nmax,
	      dissolvetime, dissolvehead, dissolvetail,
	      disstime, disshead, disstail,
	      difftime, diffhead, difftail,
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

  diffnetworktime[0] = (int)(diffhead[0]);
  diffnetworkhead[0] = (int)(diffhead[0]);
  diffnetworktail[0] = (int)(diffhead[0]);
//Rprintf("numdissolved %d numdissolve %d\n", *numdissolved, numdissolve);
  for (i = 1; i <= diffnetworkhead[0]; i++){
    diffnetworktime[i] = (int)(difftime[i]);
    diffnetworkhead[i] = (int)(diffhead[i]);
    diffnetworktail[i] = (int)(difftail[i]);
//Rprintf(" %d %d\n", dissolvedhead[i], dissolvedtail[i]);
  }

  dissnetworktime[0] = (int)(disshead[0]);
  dissnetworkhead[0] = (int)(disshead[0]);
  dissnetworktail[0] = (int)(disshead[0]);
//Rprintf("numdissolved %d numdissolve %d\n", *numdissolved, numdissolve);
  for (i = 1; i <= dissnetworkhead[0]; i++){
    dissnetworktime[i] = (int)(disstime[i]);
    dissnetworkhead[i] = (int)(disshead[i]);
    dissnetworktail[i] = (int)(disstail[i]);
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
 void MCMCSampleDynPhase2

 Using the parameters contained in the array theta, obtain the
 network statistics for a sample of size samplesize.  burnin is the
 initial number of Markov chain steps before sampling anything
 and interval is the number of MC steps between successive 
 networks in the sample.  Put all the sampled statistics into
 the networkstatistics array. 
*********************/
void MCMCSampleDynPhase2 (DynamOrder order, char *MHproposaltype, char *MHproposalpackage,
  double *theta, double gain, double *meanstats, int nphase1, int nsubphases, double *networkstatistics, 
  long int samplesize, long int burnin, 
  long int interval, int hammingterm, int fVerbose,
  double *gamma, int dyninterval,
  Edge *nmax,
  Vertex *dissolvetime, Vertex *dissolvehead, Vertex *dissolvetail,
  Vertex *disstime, Vertex *disshead, Vertex *disstail,
  Vertex *difftime, Vertex *diffhead, Vertex *difftail,
  Network *nwp, Model *m, Model *mdyn, DegreeBound *bd) {
  long int staken, tottaken, ptottaken;
  int i, j, components, diam;
  ModelTerm *mtp;
  char *fn, *sn;
  MHproposal MH;
  
//Rprintf("nsubphases %d\n", nsubphases);

  components = diam = 0;
  nwp->duration_info.MCMCtimer=0;
  
  if (fVerbose)
    Rprintf("Total m->n_stats is %i; total samplesize is %d\n",
             m->n_stats,samplesize);

  for (i = 0; MHproposaltype[i] != ' ' && MHproposaltype[i] != 0; i++);
  MHproposaltype[i] = 0;
  /* Extract the required string information from the relevant sources */
  if((fn=(char *)malloc(sizeof(char)*(i+4)))==NULL){
    error("Error in MCMCSample: Can't allocate %d bytes for fn. Memory has not been deallocated, so restart R sometime soon.\n",
	    sizeof(char)*(i+4));
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
    error("Error in ModelInitialize: Can't allocate %d bytes for sn. Memory has not been deallocated, so restart R sometime soon.\n",
	  sizeof(char)*(j+1));
  }
  sn=strncpy(sn,MHproposalpackage,i);
  sn[i]='\0';
  if (fVerbose) 
    Rprintf("MH proposal function is %s from %s package\n",fn,sn);

  /* Search for the MH proposal function pointer */
  MH.func=(void (*)(MHproposal*, DegreeBound*, Network*)) R_FindSymbol(fn,sn,NULL);
  if(MH.func==NULL){
    error("Error in MCMCSample: could not find function %s in "
	  "namespace for package %s."
	  "Memory has not been deallocated, so restart R sometime soon.\n",fn,sn);
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
  double *ubar, *u2bar, *aDdiaginv;
  ubar = (double *)malloc( m->n_stats * sizeof(double));
  u2bar = (double *)malloc( m->n_stats * sizeof(double));
  aDdiaginv = (double *)malloc( m->n_stats * sizeof(double));
  for (j=0; j < m->n_stats; j++){
    networkstatistics[j] = -meanstats[j];
    ubar[j] = 0.0;
    u2bar[j] = 0.0;
  }
  mtp = m->termarray;

  /*********************
   Burn in step.  While we're at it, use burnin statistics to 
   prepare covariance matrix for Mahalanobis distance calculations 
   in subsequent calls to M-H
   *********************/
//Rprintf("MCMCSampleDyn pre burnin numdissolve %d\n", *numdissolve);
  
    staken = 0;
    Rprintf("Starting burnin of %d steps\n", burnin);
    MetropolisHastingsDyn (order, &MH, theta,
		  networkstatistics, burnin, &staken,
		  hammingterm, fVerbose, gamma, dyninterval,
		  nmax,
		  dissolvetime, dissolvehead, dissolvetail,
		  disstime, disshead, disstail,
		  difftime, diffhead, difftail,
		  nwp, m, mdyn, bd);
    Rprintf("Phase 1: %d steps (interval = %d)\n", nphase1,burnin);
    /* Now sample networks */
    for (i=0; i <= nphase1; i++){
      MetropolisHastingsDyn (order, &MH, theta,
		  networkstatistics, burnin, &staken,
		  hammingterm, fVerbose, gamma, dyninterval,
		  nmax,
		  dissolvetime, dissolvehead, dissolvetail,
		  disstime, disshead, disstail,
		  difftime, diffhead, difftail,
		  nwp, m, mdyn, bd);
      if(i > 0){
       for (j=0; j<m->n_stats; j++){
        ubar[j]  += networkstatistics[j];
        u2bar[j] += networkstatistics[j]*networkstatistics[j];
       }
      }
//  Rprintf("done %d step gain %f \n", i, gain);
    }
    if (fVerbose){
      Rprintf("Returned from Phase 1\n");
    }
    Rprintf("\n gain times inverse variances:\n");
    for (j=0; j<m->n_stats; j++){
      aDdiaginv[j] = nphase1*gain/(u2bar[j]-ubar[j]*ubar[j]/(1.0*nphase1));
      Rprintf(" %f", aDdiaginv[j]);
    }
    Rprintf("\n");
  
    staken = 0;
    tottaken = 0;
    ptottaken = 0;
    
    /* Now sample networks */
    for (i=1; i < samplesize; i++){
      
      MetropolisHastingsDyn (order, &MH, theta,
		  networkstatistics, interval, &staken,
		  hammingterm, fVerbose, gamma, dyninterval,
		  nmax,
		  dissolvetime, dissolvehead, dissolvetail,
		  disstime, disshead, disstail,
		  difftime, diffhead, difftail,
		  nwp, m, mdyn, bd);
    /* Update theta0 */
//Rprintf("initial:\n");
      for (j=0; j<m->n_stats; j++){
        theta[j] -= aDdiaginv[j] * networkstatistics[j];
      }
//Rprintf("\n");
//    if (fVerbose){ Rprintf("nsubphases %d i %d\n", nsubphases, i); }
      if (i==(nsubphases)){
	nsubphases = trunc(nsubphases*2.52) + 1;
        if (fVerbose){Rprintf("Updating nsub to be %d\n",nsubphases);}
        for (j=0; j<m->n_stats; j++){
          aDdiaginv[j] /= 2.0;
          if (fVerbose){Rprintf("j %d theta %f ns %f\n",
		                 j, theta[j], networkstatistics[j]);}
//        if (fVerbose){ Rprintf(" %f statsmean %f",  theta[j],(networkstatistics[j]-meanstats[j])); }
        }
        Rprintf("\n");
      }
      /* Set current vector of stats equal to previous vector */
      for (j=0; j<m->n_stats; j++){
//      networkstatistics[j] -= meanstats[j];
        networkstatistics[j+m->n_stats] = networkstatistics[j];
      }
      networkstatistics += m->n_stats;
//      if (fVerbose){ Rprintf("step %d from %d:\n",i, samplesize);}
      /* This then adds the change statistics to these values */
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
  free(ubar);
  free(u2bar);
  free(MH.togglehead);
  free(MH.toggletail);
}
