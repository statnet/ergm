/*  File src/wtCD.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2017 Statnet Commons
 */
#include "wtCD.h"

/*****************
 Note on undirected networks:  For j<k, edge {j,k} should be stored
 as (j,k) rather than (k,j).  In other words, only directed networks
 should have (k,j) with k>j.
*****************/

/*****************
 void WtCD_wrapper

 Wrapper for a call from R.

 and don't forget that tail -> head
*****************/
void WtCD_wrapper(int *dnumnets, int *nedges,
		    int *tails, int *heads, double *weights, 
		    int *dn, int *dflag, int *bipartite, 
		    int *nterms, char **funnames,
		    char **sonames, 
		    char **MHproposaltype, char **MHproposalpackage,
		  double *inputs, double *theta0, int *samplesize, int *CDparams,
		  double *sample,
		    int *fVerbose, 
		    int *status){
  int directed_flag;
  Vertex n_nodes, bip, *undotail, *undohead;
  double *undoweight;
  /* Edge n_networks; */
  WtNetwork nw[1];
  WtModel *m;
  WtMHproposal MH;
  
  n_nodes = (Vertex)*dn; 
  /* n_networks = (Edge)*dnumnets;  */
  bip = (Vertex)*bipartite; 
  
  GetRNGstate();  /* R function enabling uniform RNG */
  
  directed_flag = *dflag;

  m=WtModelInitialize(*funnames, *sonames, &inputs, *nterms);

  /* Form the network */
  nw[0]=WtNetworkInitialize(tails, heads, weights, nedges[0], 
			    n_nodes, directed_flag, bip, 0, 0, NULL);

  WtMH_init(&MH,
	    *MHproposaltype, *MHproposalpackage,
	    inputs,
	    *fVerbose,
	    nw);

  undotail = calloc(MH.ntoggles * CDparams[0] * CDparams[1], sizeof(Vertex));
  undohead = calloc(MH.ntoggles * CDparams[0] * CDparams[1], sizeof(Vertex));
  undoweight = calloc(MH.ntoggles * CDparams[0] * CDparams[1], sizeof(double));
  double *extraworkspace = calloc(m->n_stats, sizeof(double));

  *status = WtCDSample(&MH,
		       theta0, sample, *samplesize, CDparams, undotail, undohead, undoweight,
		       *fVerbose, nw, m, extraworkspace);
  
  free(undotail);
  free(undohead);
  free(undoweight);
  free(extraworkspace);
  WtMH_free(&MH);
        
/* Rprintf("Back! %d %d\n",nw[0].nedges, nmax); */
  
  WtModelDestroy(m);
  WtNetworkDestroy(nw);
  PutRNGstate();  /* Disable RNG before returning */
}


/*********************
 void WtCDSample

 Using the parameters contained in the array theta, obtain the
 network statistics for a sample of size samplesize.  burnin is the
 initial number of Markov chain steps before sampling anything
 and interval is the number of MC steps between successive 
 networks in the sample.  Put all the sampled statistics into
 the networkstatistics array. 
*********************/
WtMCMCStatus WtCDSample(WtMHproposal *MHp,
			  double *theta, double *networkstatistics, 
			int samplesize, int *CDparams, Vertex *undotail, Vertex *undohead, double *undoweight, int fVerbose,
			  WtNetwork *nwp, WtModel *m, double *extraworkspace){
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

  int staken=0;
  
  /* Now sample networks */
  unsigned int i=0, sattempted=0;
  while(i<samplesize){
    
    if(WtCDStep(MHp, theta, networkstatistics, CDparams, &staken, undotail, undohead, undoweight,
		fVerbose, nwp, m, extraworkspace)!=WtMCMC_OK)
      return WtMCMC_MH_FAILED;
    
#ifdef Win32
    if( ((100*i) % samplesize)==0 && samplesize > 500){
      R_FlushConsole();
      R_ProcessEvents();
    }
#endif
      networkstatistics += m->n_stats;
      i++;

    sattempted++;
  }

  if (fVerbose){
    Rprintf("Sampler accepted %7.3f%% of %d proposed steps.\n",
	    staken*100.0/(1.0*sattempted*CDparams[0]), sattempted*CDparams[0]); 
  }
  
  return WtMCMC_OK;
}

/*********************
 void MetropolisHastings

 In this function, theta is a m->n_stats-vector just as in WtCDSample,
 but now networkstatistics is merely another m->n_stats-vector because
 this function merely iterates nsteps=CDparams[0] times through the Markov
 chain, keeping track of the cumulative change statistics along
 the way, then returns, leaving the updated change statistics in
 the networkstatistics vector.  In other words, this function 
 essentially generates a sample of size one
*********************/
WtMCMCStatus WtCDStep (WtMHproposal *MHp,
		       double *theta, double *networkstatistics,
		       int *CDparams, int *staken, Vertex *undotail, Vertex *undohead, double *undoweight,
		       int fVerbose,
		       WtNetwork *nwp,
		       WtModel *m, double *extraworkspace) {
  
  unsigned int unsuccessful=0, ntoggled=0;

  for(unsigned int step=0; step<CDparams[0]; step++){
    unsigned int mtoggled=0;
    memset(extraworkspace, 0, m->n_stats*sizeof(double));
    double cumlr = 0;
    
    for(unsigned int mult=0; mult<CDparams[1]; mult++){
      MHp->logratio = 0;
      (*(MHp->func))(MHp, nwp); /* Call MH function to propose toggles */

      if(MHp->toggletail[0]==MH_FAILED){
	switch(MHp->togglehead[0]){
	case MH_UNRECOVERABLE:
	  error("Something very bad happened during proposal. Memory has not been deallocated, so restart R soon.");
	  
	case MH_IMPOSSIBLE:
	  Rprintf("MH Proposal function encountered a configuration from which no toggle(s) can be proposed.\n");
	  return WtMCMC_MH_FAILED;
	  
	case MH_UNSUCCESSFUL:
	  warning("MH Proposal function failed to find a valid proposal.");
	  unsuccessful++;
	  if(unsuccessful>*staken*MH_QUIT_UNSUCCESSFUL){
	    Rprintf("Too many MH Proposal function failures.\n");
	    return WtMCMC_MH_FAILED;
	  }
	  continue;
	  
	case MH_CONSTRAINT:
	  MHp->logratio = -INFINITY; // Force rejection of proposal.
	  break; // Do not attempt any more proposals in this multiplicity chain.
	}
      }
      
      if(fVerbose>=5){
	Rprintf("Proposal: ");
	for(unsigned int i=0; i<MHp->ntoggles; i++)
	  Rprintf("  (%d, %d) -> %f  ", MHp->toggletail[i], MHp->togglehead[i], MHp->toggleweight[i]);
	Rprintf("\n");
      }

      /* Calculate change statistics,
	 remembering that tail -> head */
      WtChangeStats(MHp->ntoggles, MHp->toggletail, MHp->togglehead, MHp->toggleweight, nwp, m);

      // Add them to the cumulative changes.
      for(unsigned int i=0; i<m->n_stats; i++)
	extraworkspace[i] += m->workspace[i];
      
      if(fVerbose>=5){
	Rprintf("Changes: (");
	for(unsigned int i=0; i<m->n_stats; i++){
	  Rprintf(" %f ", m->workspace[i]);
	}
	Rprintf(")\n");
      }

      if(mult<CDparams[1]-1){
	/* Make proposed toggles provisionally. */
	for(unsigned int i=0; i < MHp->ntoggles; i++){
	  Vertex t=MHp->toggletail[i], h=MHp->togglehead[i];
	  double w=MHp->toggleweight[i];
	  undotail[ntoggled]=t;
	  undohead[ntoggled]=h;
	  undoweight[ntoggled]=WtGetEdge(MHp->toggletail[i], MHp->togglehead[i], nwp);
	  ntoggled++;
	  mtoggled++;
	  if(MHp->discord)
	    for(WtNetwork **nwd=MHp->discord; *nwd!=NULL; nwd++){
	      // This could be speeded up by implementing an "incrementation" function.
	      WtSetEdge(t, h, WtGetEdge(t, h, *nwd) + w - WtGetEdge(t, h, nwp), *nwd);
	    }
	  
	  WtSetEdge(t, h, w, nwp);
	}
      }

      // Accumulate the log acceptance ratio.
      cumlr += MHp->logratio;
    } // mult

    
    if(fVerbose>=5){
      Rprintf("Cumulative changes: (");
      for(unsigned int i=0; i<m->n_stats; i++)
	Rprintf(" %f ", extraworkspace[i]);
      Rprintf(")\n");
    }
    
    /* Calculate inner product */
    double ip=0;
    for (unsigned int i=0; i<m->n_stats; i++){
      ip += theta[i] * extraworkspace[i];
    }
    /* The logic is to set cutoff = ip+logratio ,
       then let the MH probability equal min{exp(cutoff), 1.0}.
       But we'll do it in log space instead.  */
    double cutoff = ip + cumlr;

    if(fVerbose>=5){
      Rprintf("log acceptance probability: %f + %f = %f\n", ip, cumlr, cutoff);
    }
    
    /* if we accept the proposed network */
    if (cutoff >= 0.0 || log(unif_rand()) < cutoff) { 
      if(fVerbose>=5){
	Rprintf("Accepted.\n");
      }
      (*staken)++; 

      if(step<CDparams[0]-1){
	/* Make the remaining proposed toggles (which we did not make provisionally) */
	for(unsigned int i=0; i < MHp->ntoggles; i++){
	  Vertex t=MHp->toggletail[i], h=MHp->togglehead[i];
	  double w=MHp->toggleweight[i];
	  undotail[ntoggled]=t;
	  undohead[ntoggled]=h;
	  undoweight[ntoggled]=WtGetEdge(MHp->toggletail[i], MHp->togglehead[i], nwp);
	  ntoggled++;

	  if(MHp->discord)
	    for(WtNetwork **nwd=MHp->discord; *nwd!=NULL; nwd++){
	      // This could be speeded up by implementing an "incrementation" function.
	      WtSetEdge(t, h, WtGetEdge(t, h, *nwd) + w - WtGetEdge(t, h, nwp), *nwd);
	    }

	  WtSetEdge(t, h, w, nwp);
	}
      }

      /* record network statistics for posterity */
      for (unsigned int i = 0; i < m->n_stats; i++){
	networkstatistics[i] += extraworkspace[i];
      }

    }else{
      if(fVerbose>=5){
	Rprintf("Rejected.\n");
      }
      // Undo the provisional toggles (the last mtoggled ones)
      for(unsigned int i=0; i < mtoggled; i++){
	ntoggled--;
	Vertex t = undotail[ntoggled], h = undohead[ntoggled];
	double w = undoweight[ntoggled];

	if(MHp->discord)
	  for(WtNetwork **nwd=MHp->discord; *nwd!=NULL; nwd++){
	    // This could be speeded up by implementing an "incrementation" function.
	    WtSetEdge(t, h, WtGetEdge(t, h, *nwd) + w - WtGetEdge(t, h, nwp), *nwd);
	  }

	WtSetEdge(t, h, w, nwp);
      }
    }
  } // step
  
  /* Undo toggles. */
  for(unsigned int i=0; i < ntoggled; i++){
    Vertex t = undotail[i], h = undohead[i];
    double w = undoweight[i];

    if(MHp->discord)
      for(WtNetwork **nwd=MHp->discord; *nwd!=NULL; nwd++){
	// This could be speeded up by implementing an "incrementation" function.
	WtSetEdge(t, h, WtGetEdge(t, h, *nwd) + w - WtGetEdge(t, h, nwp), *nwd);
      }
    
    WtSetEdge(t, h, w, nwp);
  }
  
  return WtMCMC_OK;
}

