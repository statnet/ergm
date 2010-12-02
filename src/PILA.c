#include "PILA.h"
#include "MCMC.h"
#include "R_ext/Lapack.h"
#include "R_ext/BLAS.h"

/*****************
 Note on undirected networks:  For j<k, edge {j,k} should be stored
 as (j,k) rather than (k,j).  In other words, only directed networks
 should have (k,j) with k>j.
*****************/

/*****************
 void PILA_wrapper

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
		   double *alpha, double *gamma,
		  int *fVerbose, 
                   int *attribs, int *maxout, int *maxin, int *minout,
                   int *minin, int *condAllDegExact, int *attriblength, 
		  int *mheads, int *mtails, int *mdnedges,
		  double *theta_mean_save, double *XtX_save, double *XtY_save, double *beta_save, 
		  double *direction_save, double *dtheta_save,
		  int *insensitive_save, int *ineffectual_save, int *dropped_save) {
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
	      *alpha,*gamma,
	      hammingterm,
	      (int)*fVerbose, nw, m, bd,
	      theta_mean_save, XtX_save, XtY_save, beta_save, 
	      direction_save, dtheta_save,
	      insensitive_save, ineffectual_save, dropped_save);

/*   int ii;
   double mos=0.0;
   for(ii=0; ii < bd->attrcount; ii++) 
     mos += bd->maxout[ii];
   Rprintf("bd -> attrcount = %d, sum = %f\n", ii, mos); */
        
        
/* Rprintf("Back! %d %d\n",nw[0].nedges, nmax); */
  
  ModelDestroy(m);
  if(bd)DegreeBoundDestroy(bd);
  NetworkDestroy(nw);
  if (n_medges>0 || hammingterm > 0  || formationterm > 0)
    NetworkDestroy(&nw[1]);
  PutRNGstate();  /* Disable RNG before returning */
}

double vnorm(double *a, unsigned int n){
  double tmp=0;
  for(unsigned int i=0;i<n;i++){
    tmp+=a[i]*a[i];
  }
  return(sqrt(tmp));
}

void dmat_cm_t(double *A, unsigned int n, unsigned int m, void *workspace){
  memcpy(workspace,A,n*m*sizeof(double));
  for(unsigned int r=0; r<n; r++){
    for(unsigned int c=0; c<m; c++){
      A[c+r*m]=((double *)workspace)[r+c*n];
    }
  }
}

void dmat_sq_cm_t(double *A, unsigned int n){
  for(unsigned int r=0; r<n; r++){
    for(unsigned int c=0; c<n; c++){
      double tmp=A[r+c*n];
      A[r+c*n]=A[c+r*n];
      A[c+r*n]=tmp;
    }
  }
}

void dmat_cm_delcol(double *A, unsigned int n, unsigned int m, unsigned int c){
  memmove(A+c*n,A+(c+1)*n,(m-c-1)*n*sizeof(double));
}

void dmat_cm_delrow(double *A, unsigned int n, unsigned int m, unsigned int r){
  for(unsigned int c=0; c<m; c++){
    memmove(A+c*(n-1),A+c*n,r*sizeof(double));
    memmove(A+c*(n-1)+r,A+c*n+r+1,(n-r-1)*sizeof(double));
  }
}

void PILA_SummUp(unsigned int p, unsigned int p_ext, double gamma,
		 double *theta_ext, double *cs,

		 double *theta_mean, double *n, double *cs_mag, 
		 double *XtX, double *XtY, double *cs_mean){

  (*n)*=gamma;
  (*n)++;

   for(unsigned int i=1; i<p_ext; i++){
     theta_mean[i]*=gamma;
     theta_mean[i]+=theta_ext[i]; 
   }

  for(unsigned int c=0;c<p_ext;c++){
    for(unsigned int r=0;r<p_ext;r++){
      XtX[r+c*p_ext]*=gamma;
      XtX[r+c*p_ext]+=theta_ext[r]*theta_ext[c];
    }
  }
  
  for(unsigned int c=0;c<p;c++){
    for(unsigned int r=0;r<p_ext;r++){
      XtY[r+c*p_ext]*=gamma;
      XtY[r+c*p_ext]+=theta_ext[r]*cs[c];
    }
    cs_mean[c]*=gamma;
    cs_mean[c]+=cs[c];
  }
  
  (*cs_mag)*=gamma;
  (*cs_mag)+=vnorm(cs,p);
}


#define dsaveif(var_save,var,len) {if(s%interval==0 && var_save && var) {memcpy((var_save),(var),(len)*sizeof(double)); (var_save)+=(len);}}
#define isaveif(var_save,var,len) {if(s%interval==0 && var_save && var) {memcpy((var_save),(var),(len)*sizeof(int)); (var_save)+=(len);}}
		   
/*********************
 void PILASample

 Run the joint optimization process.
*********************/
void PILASample (char *MHproposaltype, char *MHproposalpackage,
		  double *theta, double *networkstatistics, 
		 long int samplesize, unsigned int interval,
		 long int burnin,
		 double alpha, double gamma,
		  int hammingterm, int fVerbose,
		 Network *nwp, Model *m, DegreeBound *bd,
		 double *theta_mean_save, double *XtX_save, double *XtY_save, double *beta_save, 
		 double *direction_save, double *dtheta_save,
		 int *insensitive_save, int *ineffectual_save, int *dropped_save) {
  long int staken;
  MHproposal MH;
  unsigned int p_ext=m->n_stats+1, p=m->n_stats;
  
  // Init theta_ext.
  double *theta_ext=(double *)R_alloc(p_ext,sizeof(double));
  theta_ext[0]=1;
  memcpy(theta_ext+1,theta,p*sizeof(double));
  
  MH_init(&MH,
	  MHproposaltype, MHproposalpackage,
	  fVerbose,
	  nwp, bd);
  
  // Init XtX, XtY, etc.
  double *XtX = (double *)R_alloc(p_ext*p_ext,sizeof(double)),
    *XtX_work = (double *)R_alloc(p_ext*p_ext,sizeof(double)),
    *XtY = (double *)R_alloc(p*p_ext,sizeof(double)),
    *XtY_work = (double *)R_alloc(p*p_ext,sizeof(double)),
    *cs_mean = (double *)R_alloc(p,sizeof(double)),
    *ns = (double *)R_alloc(p,sizeof(double)),
    n = 0,
    *cs = (double *)R_alloc(p,sizeof(double)),
    cs_mag = 0,
    *theta_mean = (double *)R_alloc(p_ext,sizeof(double)) 
    ;

  memset(XtX,0,p_ext*p_ext*sizeof(double));
  memset(XtY,0,p*p_ext*sizeof(double));
  memset(cs_mean,0,p*sizeof(double));
  memset(cs,0,p*sizeof(double));
  memset(theta_mean,0,p_ext*sizeof(double));
  
  for(unsigned int s = 0; s<burnin; s++){
   
    memset(cs,0,p*sizeof(double));
    staken=0;
    MetropolisHastings(&MH, theta_ext+1, cs, 1, &staken,
		       hammingterm, fVerbose, nwp, m, bd);
    // Update network statistics.
    for(unsigned int i=0; i<p; i++)
      networkstatistics[i]+=cs[i];

    if(1||staken){
      PILA_SummUp(p,p_ext,gamma,
		  theta_ext,cs,

		  theta_mean,&n,&cs_mag,XtX,XtY,cs_mean);

      for(unsigned int i=0; i<p; i++)
	theta_ext[i+1]=theta[i]+rnorm(0,alpha);
    }
  }

  memcpy(theta_ext+1,theta,p*sizeof(double));
  
  unsigned int *insensitive = (unsigned int *) R_alloc(p,sizeof(unsigned int)),
    *ineffectual=(unsigned int *) R_alloc(p,sizeof(unsigned int));

  for (unsigned int s=0; s < samplesize*interval; s++){
    // Record theta 
    if(s % interval == 0){
      memcpy(theta,theta_mean+1,p*sizeof(double));
      for(unsigned int i=0; i<p;i++) theta[i]/=n;
      theta+=p;
    }
    
    memset(cs,0,p*sizeof(double));
    staken=0;
    MetropolisHastings(&MH, theta_ext+1, cs, 1, &staken,
		       hammingterm, fVerbose, nwp, m, bd);
    // Update network statistics.
    for(unsigned int i=0; i<p; i++)
      networkstatistics[i]+=cs[i];

    // Record network statistics.
    if(s && s % interval == 0){
      memcpy(networkstatistics+p,networkstatistics,p*sizeof(double));
      networkstatistics+=p;
    }
        
    if(1||staken){
      PILA_SummUp(p,p_ext,gamma,
		  theta_ext,cs,

		  theta_mean,&n,&cs_mag,XtX,XtY,cs_mean);

      // Solve for beta+:
      memcpy(XtX_work,XtX,p_ext*p_ext*sizeof(double));
      memcpy(XtY_work,XtY,p_ext*p*sizeof(double));    
      
      // "Center" the regression about the current theta_mean.
      for(unsigned int c=0;c<p_ext;c++){
	for(unsigned int r=0;r<p_ext;r++){
	  XtX_work[r+c*p_ext]-=theta_mean[r]*theta_mean[c]/n;
	}
      }
      for(unsigned int c=0;c<p;c++){
	for(unsigned int r=0;r<p_ext;r++){
	  XtY_work[r+c*p_ext]-=theta_mean[r]*cs_mean[c]/n;
	}
      }
      
      dsaveif(XtX_save,XtX,p_ext*p_ext);
      dsaveif(XtY_save,XtY,p_ext*p);
      dsaveif(theta_mean_save,theta_mean,p_ext);
    
      // Make Fortran call to solve for beta+.
      int info, p_ext_int=p_ext, p_int=p; // Last two solve annoying compiler warning
      F77_CALL(dposv)("U",&p_ext_int,&(p_int),XtX_work,
		      &p_ext_int,XtY_work,&p_ext_int,&info);
      if(info!=0) {
	Rprintf("Error in dposv, code %d, at iteration %d.",info,s);
	return;
      }
      // Probably unnecessary, but better safe than sorry:
      p_ext=p_ext_int; p=p_int;
      
      // XtY_work now contains beta+.

      dsaveif(beta_save,XtY_work,p_ext*p);
      
      // Compute the direction in which we need to go, with the approx.
      // magnitude of one step.

      memcpy(ns,networkstatistics,p*sizeof(double));
      double ns_mul=-fmin(1,cs_mag/n/vnorm(ns,p));
      for(unsigned int c=0; c<p;c++){
	ns[c]*=ns_mul;
      }
      
      dsaveif(direction_save,ns,p);
      
      // beta0 (col 1 of XtY) represent "drift" in the statistics, so
      // subtract it from the desired direction.
      for(unsigned int c=0; c<p;c++)
	ns[c]-=XtY_work[c*p_ext];

      // If a particular statistc's (dG/dtheta) is too small, 
      // mark it as being "insensitive" to
      // theta (at least given this network configuration),
      // and "drop" it.
      
      double e_max=0;
      unsigned int dropped=0;
      memset(insensitive,0,p*sizeof(double));
      // Find the max norm of the coefficient.
      for(unsigned int c=0; c<p; c++){
	if(vnorm(XtY_work+c*p_ext+1,p)>e_max){
	  e_max=vnorm(XtY_work+c*p_ext+1,p);
	}
      }
      // Mark any statistic less than the square root of machine-epsilon of it
      // as "insensitive".
      for(unsigned int c=0; c<p; c++){
	insensitive[c]=(vnorm(XtY_work+c*p_ext+1,p)/e_max<sqrt(sqrt(2.220446e-16))); // .Machine$double.eps^0.25
	if(insensitive[c]) dropped++;
      }

      isaveif(insensitive_save,insensitive,p);
      isaveif(dropped_save,&dropped,1);

      // Remove the "insensitive" statistics.
      // Counting backwards here, so need signed int.
      for(int c=p-1,deleted=0; c>=0;c--){
	if(insensitive[c]){
	  dmat_cm_delcol(XtY_work,p_ext,p-deleted,c);
	  deleted++;
	  memmove(ns+c,ns+c+1,(p-c-deleted)*sizeof(double));
	}
      }
      // XtY_work is now a matrix of dimension (p_ext)*(p-dropped).
      // For debugging purposes, zero all the "vacated" entries.
      memset(XtY_work+p_ext*(p-dropped),0,p_ext*dropped*sizeof(double));

      // Since we need to transpose beta anyway, do it here...
      // Use XtX_work for temporary storage.
      dmat_cm_t(XtY_work, p_ext, p-dropped, XtX_work);
      // beta+ is now transposed, and beta+^t has dimension
      // (p-dropped)*(p_ext).

      // Throw away the first column of beta+^t:
      dmat_cm_delcol(XtY_work,p-dropped,p_ext,0);
      // XtY now contains beta^t with dropped rows.
      // For debugging purposes, zero all the "vacated" entries.
      memset(XtY_work+p*(p-dropped),0,p*(dropped+1)*sizeof(double));

      // Also, (dropped) of the columns in the beta^t matrix will have to go to keep it square:
      // we take the norm of each column, and throw away the one with the smallest norm.
      // This is a very bad algorithm.
      memset(ineffectual,0,p*sizeof(unsigned int));
      for(unsigned int i=0; i<dropped;i++){
	unsigned int c_min=0;
	double m_min=HUGE_VAL;
	for(unsigned int c=0; c<p; c++){
	  if(!ineffectual[c] && vnorm(XtY_work+c*(p-dropped),p-dropped)<m_min){
	    c_min=c;
	    m_min=vnorm(XtY_work+c*(p-dropped),p-dropped);
	  }
	}
	ineffectual[c_min]=1;
      }
      for(int c=p-1,deleted=0; c>=0;c--){
	if(ineffectual[c]){
	  dmat_cm_delcol(XtY_work,p-dropped,p-deleted,c);
	  deleted++;
	}
      }
      
      isaveif(ineffectual_save,ineffectual,p);
      
      // XtY_work is now a p-dropped by p-dropped matrix,
      // whose rows correspond to network statistics and whose columns
      // correspond to parameters.
      // For debugging purposes, zero all the "vacated" entries.
      memset(XtY_work+(p-dropped)*(p-dropped),0,(p*p_ext-(p-dropped)*(p-dropped))*sizeof(double));
    
      // Solve for G * beta^-1.
      int one=1;
      int msize = p-dropped;
      // Note that XtX_work is throw-away memory here.
      F77_CALL(dgesv)(&msize,&one,XtY_work,&msize,
		      (int *) XtX_work,ns,&msize,&info);
      if(info!=0) {
	Rprintf("Error in dgesv, code %d, at iteration %d.",info,s);
	return;
      }
      // ns now contains the estimated change in theta with "ineffectual"
      // coefficients dropped. Now, reconstruct the full-length change in theta.
      double stepsize=vnorm(ns,p-dropped)/(p-dropped);
      for(int i=p-1,pos=p-dropped-1; i>=0; i--){
	if(ineffectual[i])
	  ns[i]=rnorm(0,stepsize/3);
	else{
	  // Note that pos <= i, always.
	  ns[i]=ns[pos];
	  pos--;
	}
      }
      
      // Now, generate new theta by adding the change to the mean.
      
      // Normalize the change to have magnitude alpha:
      
      double stepscl=alpha/vnorm(ns,p);

      dsaveif(dtheta_save,ns,p);

      for(unsigned int i=0; i<p;i++)
	theta_ext[i+1]=theta_mean[i+1]/n+ns[i]*stepscl;
    }
  }
  MH_free(&MH);
}

#undef dsaveif
#undef isaveif
