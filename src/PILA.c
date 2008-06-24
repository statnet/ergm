#include "PILA.h"
#include "R_ext/Lapack.h"
#include "R_ext/BLAS.h"

/*****************
 Note on undirected networks:  For j<k, edge {j,k} should be stored
 as (j,k) rather than (k,j).  In other words, only directed networks
 should have (k,j) with k>j.
*****************/

/*****************
 void MCMC_wrapper

 Wrapper for a call from R.
*****************/
void PILA_wrapper(int *heads, int *tails, int *dnedges,
		  int *dn, int *dflag, int *bipartite, 
                   int *nterms, char **funnames,
                   char **sonames, 
                   char **MHproposaltype, char **MHproposalpackage,
                   double *inputs, double *theta,
		   int *samplesize, int *interval, 
                   double *sample, int *burnin, 
		   double *theta_burnin, double *sample_burnin,
		   double *alpha, double *gamma,
		  int *fVerbose, 
                   int *attribs, int *maxout, int *maxin, int *minout,
                   int *minin, int *condAllDegExact, int *attriblength, 
		  int *mheads, int *mtails, int *mdnedges) {
  int directed_flag, hammingterm, formationterm;
  Vertex n_nodes, nmax, bip, hhead, htail;
  Edge n_edges, n_medges, nddyads, kedge;
  Network nw[2];
  DegreeBound *bd;
  Model *m;
  ModelTerm *thisterm;
  
  n_nodes = (Vertex)*dn; /* coerce double *dn to type Vertex */
  n_edges = (Edge)*dnedges; /* coerce int *dnedges to type Edge */
  n_medges = 0; //(Edge)*mdnedges; /* coerce int *mdnedges to type Edge --- turned off for now */
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
  formationterm=ModelTermFormation (*funnames, *nterms);
  if(formationterm>0){
   Network nwformation;
   thisterm = m->termarray + formationterm - 1;
   nddyads = (Edge)(thisterm->inputparams[0]);
   nwformation=NetworkInitializeD(thisterm->inputparams+1,
				  thisterm->inputparams+1+nddyads, nddyads,
          n_nodes, directed_flag, bip,0);
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
  
  bd=DegreeBoundInitialize(attribs, maxout, maxin, minout, minin,
			   *condAllDegExact, *attriblength, nw);
  PILASample (*MHproposaltype, *MHproposalpackage,
	      theta, sample, (long int)*samplesize,
	       *interval,
	      (long int)*burnin,
	       theta_burnin, sample_burnin, 
	      *alpha,*gamma,
	      hammingterm,
	      (int)*fVerbose, nw, m, bd);

/*   int ii;
   double mos=0.0;
   for(ii=0; ii < bd->attrcount; ii++) 
     mos += bd->maxout[ii];
   Rprintf("bd -> attrcount = %d, sum = %f\n", ii, mos); */
        
        
/* Rprintf("Back! %d %d\n",nw[0].nedges, nmax); */
  
  ModelDestroy(m);
  DegreeBoundDestroy(bd);
  NetworkDestroy(nw);
  if (n_medges>0 || hammingterm > 0  || formationterm > 0)
    NetworkDestroy(&nw[1]);
  PutRNGstate();  /* Disable RNG before returning */
}

double pythagn(double *a, unsigned int n){
  double tmp=0;
  for(unsigned int i=0;i<n;i++){
    tmp+=a[i]*a[i];
  }
  return(sqrt(tmp));
}

/*********************
 void PILASample

 Run the joint optimization process.
*********************/
void PILASample (char *MHproposaltype, char *MHproposalpackage,
		  double *theta, double *networkstatistics, 
		  long int samplesize, unsigned int interval,
		  long int burnin,
		  double *theta_burnin, double *sample_burnin,
		  double alpha, double gamma,
		  int hammingterm, int fVerbose,
		  Network *nwp, Model *m, DegreeBound *bd) {
  long int staken;
  MHproposal MH;
  unsigned int p_ext=m->n_stats+1;
  
  // Init theta_ext.
  double *theta_ext=(double *)R_alloc(p_ext,sizeof(double));
  theta_ext[0]=1;
  memcpy(theta_ext+1,theta,m->n_stats*sizeof(double));
  
  MH_init(&MH,
	  MHproposaltype, MHproposalpackage,
	  fVerbose,
	  nwp, bd);
  
  // Init XtX, XtY, etc.
  double *XtX = (double *)R_alloc(p_ext*p_ext,sizeof(double)),
    *XtX_work = (double *)R_alloc(p_ext*p_ext,sizeof(double)),
    *XtY = (double *)R_alloc(m->n_stats*p_ext,sizeof(double)),
    *XtY_work = (double *)R_alloc(m->n_stats*p_ext,sizeof(double)),
    *ns_dev = (double *)R_alloc(m->n_stats,sizeof(double)),
    *ns_dev_work = (double *)R_alloc(m->n_stats,sizeof(double)),
    n = 0,
    *cs = (double *)R_alloc(m->n_stats,sizeof(double)),
    cs_mag = 0,
    *theta_mean = (double *)R_alloc(p_ext,sizeof(double)) 
    ;
  memset(XtX,0,p_ext*p_ext*sizeof(double));
  memset(XtY,0,m->n_stats*p_ext*sizeof(double));
  memset(ns_dev,0,m->n_stats*sizeof(double));
  memset(cs,0,m->n_stats*sizeof(double));
  
  for(unsigned int s = 0; s<burnin; s++){
    n*=gamma;
    n++;

    // Record theta for burnin
    if(theta_burnin){
      memcpy(theta_burnin,theta_ext+1,m->n_stats*sizeof(double));
      theta_burnin+=m->n_stats;
    }

    // Update mean theta (the 1st element stays at 0)
    for(unsigned int i=1; i<p_ext; i++){
      theta_mean[i]*=gamma;
      theta_mean[i]+=theta_ext[i]; 
    }

    memset(cs,0,m->n_stats*sizeof(double));
    MetropolisHastings(&MH, theta_ext+1, cs, 1, &staken,
		       hammingterm, fVerbose, nwp, m, bd);
    // Update network statistics.
    for(unsigned int i=0; i<m->n_stats; i++)
      networkstatistics[i]+=cs[i];
    // Record network statistics for burnin.
    if(sample_burnin){
      memcpy(sample_burnin,networkstatistics,m->n_stats*sizeof(double));
      sample_burnin+=m->n_stats;
    }
    // update XtX, XtY, and ns_dev
    for(unsigned int c=0;c<p_ext;c++){
      for(unsigned int r=0;r<p_ext;r++){
	XtX[r+c*p_ext]*=gamma;
	XtX[r+c*p_ext]+=(theta_ext[r]-theta_mean[r]/n)*(theta_ext[c]-theta_mean[c]/n);
      }
    }
    
    for(unsigned int c=0;c<m->n_stats;c++){
      for(unsigned int r=0;r<p_ext;r++){
	XtY[r+c*p_ext]*=gamma;
	XtY[r+c*p_ext]+=(theta_ext[r]-theta_mean[r]/n)*cs[c];
      }
      ns_dev[c]*=gamma;
      ns_dev[c]+=networkstatistics[c];
    }
    
    cs_mag*=gamma;
    cs_mag+=pythagn(cs,m->n_stats);
    
    for(unsigned int i=1; i<p_ext; i++)
      theta_ext[i]=theta[i-1]+rnorm(0,.5);
  }

  memcpy(theta_ext+1,theta,m->n_stats*sizeof(double));
  for (unsigned int s=0; s < samplesize*interval; s++){
    // Record theta 
    if(s % interval == 0){
      memcpy(theta,theta_mean+1,m->n_stats*sizeof(double));
      for(unsigned int i=0; i<m->n_stats;i++) theta[i]/=n;
      theta+=m->n_stats;
    }
    
    n*=gamma;
    n++;

    // Update mean theta (the 1st element stays at 0)
    
    for(unsigned int i=1; i<p_ext; i++){
      theta_mean[i]*=gamma;
      theta_mean[i]+=theta_ext[i]; 
    }

    memset(cs,0,m->n_stats*sizeof(double));
    MetropolisHastings(&MH, theta_ext+1, cs, 1, &staken,
		       hammingterm, fVerbose, nwp, m, bd);
    // Update network statistics.
    for(unsigned int i=0; i<m->n_stats; i++)
      networkstatistics[i]+=cs[i];
    // Record network statistics.
    if(s % interval == 0){
      memcpy(networkstatistics+m->n_stats,networkstatistics,m->n_stats*sizeof(double));
      networkstatistics+=m->n_stats;
    }

    // update XtX, XtY, and ns_dev
    for(unsigned int c=0;c<p_ext;c++){
      for(unsigned int r=0;r<p_ext;r++){
	XtX[r+c*p_ext]*=gamma;
	XtX[r+c*p_ext]+=(theta_ext[r]-theta_mean[r]/n)*(theta_ext[c]-theta_mean[c]/n);
      }
    }
    
    for(unsigned int c=0;c<m->n_stats;c++){
      for(unsigned int r=0;r<p_ext;r++){
	XtY[r+c*p_ext]*=gamma;
	XtY[r+c*p_ext]+=(theta_ext[r]-theta_mean[r]/n)*cs[c];
      }
      ns_dev[c]*=gamma;
      ns_dev[c]+=networkstatistics[c];
    }
    
    cs_mag*=gamma;
    cs_mag+=pythagn(cs,m->n_stats);
        
    memcpy(XtX_work,XtX,p_ext*p_ext*sizeof(double));
    memcpy(XtY_work,XtY,p_ext*m->n_stats*sizeof(double));
    memcpy(ns_dev_work,ns_dev,m->n_stats*sizeof(double));
    
    //Rprintf("\n%u:\n",s);
    //Rprintf("cs: %f\n",cs[0]);
    //Rprintf("ns: %f\n",networkstatistics[0]);
    //Rprintf("theta_mean: %f\n",theta_mean[1]/n);
    //Rprintf("XtX: %f ,%f ,%f ,%f\n",XtX[0]/n,XtX[1]/n,XtX[2]/n,XtX[3]/n);
    //Rprintf("XtY: %f ,%f\n",XtY[0]/n,XtY[1]/n);
    //Rprintf("ns_dev: %f\n",ns_dev[0]/n);
    // Solve for beta+.
    int info;
    F77_CALL(dposv)("U",&p_ext,&(m->n_stats),XtX_work,
		    &p_ext,XtY_work,&p_ext,&info);
    if(info!=0) {error("Error in dposv, code %d.",info);}
    // Now, XtY_work contains beta+.
    //Rprintf("beta+: %f ,%f\n",XtY_work[0],XtY_work[1]);
    // Compute the vector with direction we need to head in,
    // and having the magnitude of a one-step change.
    // Note that this is a simplification, since pythagn(ns_dev)
    // is n times the magnitude it should be, as is ns_dev_work.
    for(unsigned int c=0; c<m->n_stats;c++)
      ns_dev_work[c]*=(cs_mag/n)/pythagn(ns_dev,m->n_stats);
    //Rprintf("Target CS (w/o beta0): %f\n",ns_dev_work[0]);
    // beta0 represent "drift" in the statistics, and so should be
    // subtracted from the desired direction
    for(unsigned int c=0; c<m->n_stats;c++)
      ns_dev_work[c]+=XtY_work[c*m->n_stats];
    //Rprintf("Target CS (w/ beta0): %f\n",ns_dev_work[0]);
    // The first row of beta needs to be thrown
    // away before inverting.
    for(unsigned int c=0;c<m->n_stats;c++){
      for(unsigned int r=1;r<p_ext;r++){
	XtY_work[r+c*m->n_stats-1]=XtY_work[r+c*p_ext];
      }
    }

    // In order to write the problem as Ax=b, transpose beta.
    for(unsigned int c=0;c<m->n_stats;c++){
      for(unsigned int r=0;r<m->n_stats;r++){
	double tmp=XtY_work[r+c*m->n_stats];
	XtY_work[r+c*m->n_stats]=XtY_work[c+r*m->n_stats];
	XtY_work[c+r*m->n_stats]=tmp;
      }
    }
    
    //Rprintf("beta: %f\n",XtY_work[0]);
    
    // Solve for G * beta^-1.
    int one=1;
    // Note that XtX_work is throw-away memory here.
    F77_CALL(dgesv)(&(m->n_stats),&one,XtY_work,&(m->n_stats),
		    (int *) XtX_work,ns_dev_work,&(m->n_stats),&info);
    if(info!=0) {error("Error in dgesv, code %d.",info);}
    // ns_dev_work now contains the estimated change in theta.
    //Rprintf("change in theta (raw): %f\n",ns_dev_work[0]);
    // Now, generate new theta by adding the change to the mean,
    // along with a random perturbation.
    for(unsigned int i=0; i<m->n_stats;i++)
      theta_ext[i+1]=theta_mean[i+1]/n-alpha*ns_dev_work[i]+rnorm(0,0.5);
  }
  MH_free(&MH);
  
}

