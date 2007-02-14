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
void MCMC_wrapper (double *heads, double *tails, double *dnedges,
                   double *dn, int *dflag, double *bipartite, 
                   int *nterms, char **funnames,
                   char **sonames, 
                   char **MHproposaltype, char **MHproposalpackage,
                   double *inputs, double *theta0, double *samplesize, 
                   double *sample, double *burnin, double *interval,  
                   int *newnetwork, 
                   int *fVerbose, 
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
  Model *m;
  ModelTerm *thisterm;
  
  n_nodes = (Vertex)*dn; /* coerce double *dn to type Vertex */
  n_edges = (Vertex)*dnedges; /* coerce double *dnedges to type Vertex */
  n_medges = (Vertex)*mdnedges; /* coerce double *mdnedges to type Vertex */
  nmax = (Vertex)*maxedges; /* coerce double *maxedges to type Vertex */
  bip = (Vertex)*bipartite; /* coerce double *bipartite to type Vertex */
  
  GetRNGstate();  /* R function enabling uniform RNG */
  
  directed_flag = *dflag;

  for (i = 0; i < nmax; i++)
    newnetwork[i] = 0;

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
  formationterm=ModelTermFormation (*funnames, *nterms);
  if(formationterm>0){
   Network nwformation;
   thisterm = m->termarray + formationterm - 1;
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
  MCMCSample (*MHproposaltype, *MHproposalpackage,
	      theta0, sample, (long int)*samplesize,
	      (long int)*burnin, (long int)*interval,
	      hammingterm,
	      (int)*fVerbose, nw, m, bd);
  
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
          newnetwork[nextedge] = v;
	  nextedge++;
          newnetwork[nextedge] = nw[0].outedges[e].value;
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
      newnetwork[nextedge] = k;
      nextedge++;
      newnetwork[nextedge] = v;
      nextedge++;
	  }else{
      newnetwork[nextedge] = v;
      nextedge++;
      newnetwork[nextedge] = k;
      nextedge++;
	  }
	}
     }
  }
  newnetwork[0]=nextedge;

  ModelDestroy(m);
  DegreeBoundDestroy(bd);
  NetworkDestroy(nw);
  if (n_medges>0 || hammingterm > 0  || formationterm > 0)
    NetworkDestroy(&nw[1]);
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
void MCMCSample (char *MHproposaltype, char *MHproposalpackage,
  double *theta, double *networkstatistics, 
  long int samplesize, long int burnin, 
  long int interval, int hammingterm, int fVerbose,
  Network *nwp, Model *m, DegreeBound *bd) {
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
  MetropolisHastings(&MH, theta, networkstatistics, burnin, &staken,
		     hammingterm, fVerbose, nwp, m, bd);
  
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
      
      MetropolisHastings (&MH, theta, networkstatistics, interval, &staken,
			  hammingterm, fVerbose, nwp, m, bd);
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
    }
    /*********************
    Below is an extremely crude device for letting the user know
    when the chain doesn't accept many of the proposed steps.
    *********************/
    if (fVerbose){
      Rprintf("Metropolis-Hastings accepted %7.3f%% of %d steps.\n",
	      tottaken*100.0/(1.0*interval*samplesize), interval*samplesize); 
    }
  }else{
    if (fVerbose){
      Rprintf("Metropolis-Hastings accepted %7.3f%% of %d steps.\n",
	      staken*100.0/(1.0*burnin), burnin); 
    }
  }
  free(MH.togglehead);
  free(MH.toggletail);
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
			 long int nsteps, long int *staken,
			 int hammingterm, int fVerbose,
			 Network *nwp,
       Model *m, DegreeBound *bd) {
  long int step, taken;
  int i, curstat=0;
  double *dstats, ip, cutoff;
  ModelTerm *mtp;
  
  step = taken = 0;
/*  if (fVerbose)
    Rprintf("Now proposing %d MH steps... ", nsteps); */
  while (step < nsteps) {
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
        ToggleEdgeWithTimestamp(MHp->togglehead[i], MHp->toggletail[i], nwp);
        if(hammingterm){
	 ToggleEdge(MHp->togglehead[i],  MHp->toggletail[i], &nwp[1]);  /* Toggle the discord for this edge */
	}
      }
      /* record network statistics for posterity */
      for (i = 0; i < m->n_stats; i++)
        networkstatistics[i] += m->workspace[i];
      taken++;
    }
    step++;
    nwp->duration_info.MCMCtimer++;
  }
/*  if (fVerbose)
    Rprintf("%d taken (MCMCtimer=%d)\n", taken, nwp->duration_info.MCMCtimer); */

  *staken = taken;
}

/********************
 int CheckTogglesValid
********************/
int CheckTogglesValid(MHproposal *MHp, DegreeBound *bd, Network *nwp) {
  int fvalid;
  int i;
  int *hattr = (int *) malloc(sizeof(int) * bd->attrcount);
  int *tattr = (int *) malloc(sizeof(int) * bd->attrcount);

  fvalid = 1;

  /* Make proposed toggles */
  for (i=0; i<MHp->ntoggles; i++)
    ToggleEdge(MHp->togglehead[i], MHp->toggletail[i], nwp);

  /* if we're bounding degrees by attribute */
  if (bd->fBoundDegByAttr && fvalid)
    {      
      Edge e;
      Vertex v;
      int k;
      
      if (nwp->directed_flag)
	{
	  /* for each head and tail pair */
	  for (i = 0; i < MHp->ntoggles && fvalid; i++)
	    {
              for (k=0; k < bd->attrcount; k++){
	        hattr[k] = tattr[k] = 0;
	      }
	      /* calculate head outdegree totals for each attribute
		 for each outedge of the head 	      */
	      
	      for(e = EdgetreeMinimum(nwp->outedges, MHp->togglehead[i]);
		  (v = nwp->outedges[e].value) != 0;
		  e = EdgetreeSuccessor(nwp->outedges, e))
		{
		  for (k=0; k < bd->attrcount; k++)
		    if (bd->attribs[v-1 + k*nwp->nnodes]) hattr[k]++;
		}
	      
	      /* calculate tail indegree totals for each attribute
		 for each inedge of the tail */
	      
	      for(e = EdgetreeMinimum(nwp->inedges, MHp->toggletail[i]);
		  (v = nwp->inedges[e].value) != 0;
		  e = EdgetreeSuccessor(nwp->inedges, e))
		{
		  for (k=0; k < bd->attrcount; k++)
		    if (bd->attribs[v-1 + k*nwp->nnodes]) tattr[k]++;
		}

	      /* for each attribute */


	      for (k=0; k < bd->attrcount && fvalid; k++){
		fvalid=!((hattr[k]>bd->maxout[MHp->togglehead[i]-1+k*nwp->nnodes])||
		  (hattr[k] < bd->minout[MHp->togglehead[i]-1+k*nwp->nnodes]) || 
		  (tattr[k] >  bd->maxin[MHp->toggletail[i]-1+k*nwp->nnodes]) ||
		  (tattr[k] <  bd->minin[MHp->toggletail[i]-1+k*nwp->nnodes]) );
	      }
	    }
	}
      else /* ! nwp->directed_flag  */
	{
	  /* for each head and tail pair */
	  for (i = 0; i < MHp->ntoggles && fvalid; i++)
	    {
              for (k=0; k < bd->attrcount; k++){
	        hattr[k] = tattr[k] = 0;
	      }
	      
	      /* calculate head totals for each attribute
		 for each outedge and inedge of the head  */
	      
	      for(e = EdgetreeMinimum(nwp->outedges, MHp->togglehead[i]);
		  (v = nwp->outedges[e].value) != 0;
		  e = EdgetreeSuccessor(nwp->outedges, e))
		{
		  for (k=0; k < bd->attrcount; k++)
		    if (bd->attribs[v-1 + k*nwp->nnodes])
		      hattr[k]++;
		}
	      for(e = EdgetreeMinimum(nwp->inedges, MHp->togglehead[i]);
		  (v = nwp->inedges[e].value) != 0;
		  e = EdgetreeSuccessor(nwp->inedges, e))
		{
		  for (k=0; k < bd->attrcount; k++)
		    if (bd->attribs[v-1 + k*nwp->nnodes])
		      hattr[k]++;
		}
	      
	      /* calculate tail totals for each attribute
		 for each outedge and inedge of the tail */
	      
	      for(e = EdgetreeMinimum(nwp->outedges, MHp->toggletail[i]);
		  (v = nwp->outedges[e].value) != 0;
		  e = EdgetreeSuccessor(nwp->outedges, e))
		{
		  for (k=0; k < bd->attrcount; k++)
		    if (bd->attribs[v-1 + k*nwp->nnodes])
		      tattr[k]++;
		}
	      for(e = EdgetreeMinimum(nwp->inedges, MHp->toggletail[i]);
		  (v = nwp->inedges[e].value) != 0;
		  e = EdgetreeSuccessor(nwp->inedges, e))
		{
		  for (k=0; k < bd->attrcount; k++)
		    if (bd->attribs[v-1 + k*nwp->nnodes])
		      tattr[k]++;
		}

	      /* for each attribute
		 check heads' and tails' outmax and outmin */
	      for (k=0; k < bd->attrcount && fvalid; k++){
		fvalid=!((hattr[k]>bd->maxout[MHp->togglehead[i]-1+k*nwp->nnodes])|| 
		  (hattr[k] < bd->minout[MHp->togglehead[i]-1+k*nwp->nnodes]) || 
		  (tattr[k] > bd->maxout[MHp->toggletail[i]-1+k*nwp->nnodes]) ||
		  (tattr[k] < bd->minout[MHp->toggletail[i]-1+k*nwp->nnodes]) );
	      }
	    }
	}
    }
  
  free(hattr);
  free(tattr);

  /* Make proposed toggles */
  for (i=0; i<MHp->ntoggles; i++)
    ToggleEdge(MHp->togglehead[i], MHp->toggletail[i], nwp);

  return fvalid;
}

int CheckConstrainedTogglesValid(MHproposal *MHp, DegreeBound *bd, Network *nwp)
{
  int fvalid = 1;
  int i;

  /* Make proposed toggles */
  for (i=0; i<MHp->ntoggles; i++)
    ToggleEdge(MHp->togglehead[i], MHp->toggletail[i], nwp);

  /* if we're bounding degrees by attribute */
  if (bd->fBoundDegByAttr && fvalid)
    {
      
      Edge e;
      Vertex v;
      int k;
      int *hattr = (int *) malloc(sizeof(int) * bd->attrcount);
      int *tattr = (int *) malloc(sizeof(int) * bd->attrcount);
      
      if (nwp->directed_flag)
	{
	  /* for each head and tail pair */
	  for (i = 0; i < MHp->ntoggles && fvalid; i++)
	    {
              for (k=0; k < bd->attrcount; k++){
	        hattr[k] = tattr[k] = 0;
	      }
	      /* calculate head outdegree totals for each attribute
		 for each outedge of the head 	      */
	      
	      for(e = EdgetreeMinimum(nwp->outedges, MHp->togglehead[i]);
		  (v = nwp->outedges[e].value) != 0;
		  e = EdgetreeSuccessor(nwp->outedges, e))
		{
		  for (k=0; k < bd->attrcount; k++)
		    if (bd->attribs[v-1 + k*nwp->nnodes])
		      hattr[k]++;
		}
	      
	      /* calculate tail indegree totals for each attribute
		 for each inedge of the tail */
	      
	      for(e = EdgetreeMinimum(nwp->inedges, MHp->toggletail[i]);
		  (v = nwp->inedges[e].value) != 0;
		  e = EdgetreeSuccessor(nwp->inedges, e))
		{
		  for (k=0; k < bd->attrcount; k++)
		    if (bd->attribs[v-1 + k*nwp->nnodes])
		      tattr[k]++;
		}

	      /* for each attribute */
	      for (k=0; k < bd->attrcount && fvalid; k++){
		fvalid=!((hattr[k]>bd->maxout[MHp->togglehead[i]-1+k*nwp->nnodes])||
		  (hattr[k] < bd->minout[MHp->togglehead[i]-1+k*nwp->nnodes]) || 
		  (tattr[k] >  bd->maxin[MHp->toggletail[i]-1+k*nwp->nnodes]) ||
		  (tattr[k] <  bd->minin[MHp->toggletail[i]-1+k*nwp->nnodes])) ;
	      }
	    }
	}
      else /* ! nwp->directed_flag */
	{
	  /* for each head and tail pair */
	  for (i = 0; i < MHp->ntoggles && fvalid; i++)
	    {
              for (k=0; k < bd->attrcount; k++){
	        hattr[k] = tattr[k] = 0;
	      }
	      
	      /* calculate head totals for each attribute
		 for each outedge and inedge of the head  */
	      
	      for(e = EdgetreeMinimum(nwp->outedges, MHp->togglehead[i]);
		  (v = nwp->outedges[e].value) != 0;
		  e = EdgetreeSuccessor(nwp->outedges, e))
		{
		  for (k=0; k < bd->attrcount; k++)
		    if (bd->attribs[v-1 + k*nwp->nnodes])
		      hattr[k]++;
		}
	      for(e = EdgetreeMinimum(nwp->inedges, MHp->togglehead[i]);
		  (v = nwp->inedges[e].value) != 0;
		  e = EdgetreeSuccessor(nwp->inedges, e))
		{
		  for (k=0; k < bd->attrcount; k++)
		    if (bd->attribs[v-1 + k*nwp->nnodes])
		      hattr[k]++;
		}
	      
	      /* calculate tail totals for each attribute
		 for each outedge and inedge of the tail */
	      for(e = EdgetreeMinimum(nwp->outedges, MHp->toggletail[i]);
		  (v = nwp->outedges[e].value) != 0;
		  e = EdgetreeSuccessor(nwp->outedges, e))
		{
		  for (k=0; k < bd->attrcount; k++)
		    if (bd->attribs[v-1 + k*nwp->nnodes])
		      tattr[k]++;
		}
	      for(e = EdgetreeMinimum(nwp->inedges, MHp->toggletail[i]);
		  (v = nwp->inedges[e].value) != 0;
		  e = EdgetreeSuccessor(nwp->inedges, e))
		{
		  for (k=0; k < bd->attrcount; k++)
		    if (bd->attribs[v-1 + k*nwp->nnodes])
		      tattr[k]++;
		}

	      /* for each attribute
		 check heads' and tails' outmax and outmin */
	      for (k=0; k < bd->attrcount && fvalid; k++)
		fvalid=!(hattr[k]>bd->maxout[MHp->togglehead[i]-1+k*nwp->nnodes])||
		  (hattr[k] < bd->minout[MHp->togglehead[i]-1+k*nwp->nnodes]) || 
		  (tattr[k] > bd->maxout[MHp->toggletail[i]-1+k*nwp->nnodes]) ||
		  (tattr[k] < bd->minout[MHp->toggletail[i]-1+k*nwp->nnodes]) ;
	    }
	}
      free(hattr);
      free(tattr);
    }
    /* Make proposed toggles */
  for (i=0; i<MHp->ntoggles; i++)
    ToggleEdge(MHp->togglehead[i], MHp->toggletail[i], nwp);
  
  return fvalid;
}

/*****************
 void MCMC_global

 Wrapper for a call from R.  Return the change in the statistics when
 we go from the observed graph to an empty graph.  If the empty graph
 has true global values equal to zero for all statistics, then this change
 (with a minus sign in front) gives the true global values for the
 observed graph.
*****************/
void MCMC_global (double *heads, double *tails, double *dnedges,
		  double *dn, int *dflag,  double *bipartite,
		  int *nterms, char **funnames,
		  char **sonames, double *inputs,  double *stats)
{	
  int i, directed_flag, hammingterm, formationterm;
  Vertex n_nodes, hhead, htail;
  Edge n_edges, nddyads, kedge;
  Network nw[2];
  Model *m;
  ModelTerm *mtp;
  ModelTerm *thisterm;
  double *dstats;
  Vertex *head, *tail;
  Vertex bip;

//	     Rprintf("prestart with setup\n");
  n_nodes = (Vertex)*dn; 
  n_edges = (Edge)*dnedges;     
  directed_flag = *dflag;
  bip = (Vertex)*bipartite; /* coerce double *bipartite to type Vertex */
  
  head = (Vertex *) malloc(sizeof(Vertex) * n_edges);
  tail = (Vertex *) malloc(sizeof(Vertex) * n_edges);
  for (i = 0; i < n_edges; i++) { /* coerce edgelist to Vertex from double */
    if ( !(directed_flag) && heads[i] > tails[i] ) {
      head[i] = (Vertex)tails[i];
      tail[i] = (Vertex)heads[i];
    } else {
      head[i] = (Vertex)heads[i];
      tail[i] = (Vertex)tails[i];
    }
  }
  
  m=ModelInitialize(*funnames, *sonames, inputs, *nterms);
  nw[0]=NetworkInitialize(heads, tails, n_edges, n_nodes, directed_flag, bip);

  hammingterm=ModelTermHamming (*funnames, *nterms);
//	     Rprintf("start with setup\n");
  if(hammingterm>0){
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
   for (kedge=1; kedge <= nwhamming.nedges; kedge++) {
     FindithEdge(&hhead, &htail, kedge, &nwhamming);
     if(EdgetreeSearch(hhead, htail, nw[0].outedges) == 0){
       ToggleEdge(hhead, htail, &nw[1]);
     }
   }
   for (kedge=1; kedge <= nw[0].nedges; kedge++) {
     FindithEdge(&hhead, &htail, kedge, &nw[0]);
     if(EdgetreeSearch(hhead, htail, nwhamming.outedges) == 0){
       ToggleEdge(hhead, htail, &nw[1]);
     }
   }
   free(dhead);
   free(dtail);
   NetworkDestroy(&nwhamming);
//	     Rprintf("done with setup nw[1].nedges %d\n",nw[1].nedges);
  }

// Really this is a formation term
  formationterm=ModelTermFormation (*funnames, *nterms);
  if(formationterm>0){
   Network nwformation;
   thisterm = m->termarray + formationterm - 1;
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

  mtp = m->termarray; /* points to first model term for now */
  dstats = stats; /* points to start of change-stats vector to be returned */
  for (i=0; i < m->n_stats; i++) 
    dstats[i]=0.0; /* Initialize to zero */
  
  for (i=0; i < *nterms; i++) { /* Calculate change statistics.  Note that
    we're sending as toggles the entire list of edges.  These toggles will
    take us from the observed graph to the empty graph.  */
    mtp->dstats = dstats;
    (*(mtp->func))(n_edges, head, tail, mtp, nw);
    dstats += (mtp++)->nstats; /* increment mtp and dstats to next ModelTerm */
  }

  ModelDestroy(m);
  NetworkDestroy(nw);
  if (hammingterm > 0 || formationterm > 0)
    NetworkDestroy(&nw[1]);
  free(head);
  free(tail);
}

