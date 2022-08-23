/*  File src/MPLE.c in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2022 Statnet Commons
 */
#include "MPLE.h"
#include "ergm_changestat.h"
#include "ergm_rlebdm.h"

static double **MPLE_workspace = NULL;
static StoreDVecMapENE *MPLE_covfreq = NULL;
static khint_t MPLE_nalloc = 0;
static khint_t MPLE_nalloc_max = 0;

// TODO: Consider preallocating space for these vectors parcelling them out.
static inline double *MPLE_workspace_push(double *ptr){
  if(MPLE_nalloc == MPLE_nalloc_max)
    MPLE_workspace = Realloc(MPLE_workspace, (MPLE_nalloc_max = MAX(MPLE_nalloc_max, 1u) * 2u), double*);

  MPLE_workspace[MPLE_nalloc++] = ptr;
  return ptr;
}

SEXP MPLE_workspace_free(){
  if(MPLE_covfreq){
    kh_destroy(DVecMapENE, MPLE_covfreq);
    MPLE_covfreq = NULL;
  }

  if(MPLE_workspace){
    for(unsigned int i = 0; i < MPLE_nalloc; i++) Free(MPLE_workspace[i]);
    Free(MPLE_workspace);
    MPLE_nalloc_max = MPLE_nalloc = 0;
  }

  return R_NilValue;
}

/* *****************
 void MPLE_wrapper

 Wrapper for a call from R.  
 Only find the MPLE, so no MCMC stuff is necessary. 

 Re-Written by David Schruth to add compression:
   Since there are often many different 0->1 dyad changes that
   result in identical change statistic vectors, this routine
   now returns only the unique change statistic vectors in the 
   covmat matrix, along with a separate vector, weightsvector,
   giving the number of repetitions of each unique set of change
   statistics.  Note that two sets of statistics are considered
   unique in this context if the corresponding response values
   (i.e., dyad values, edge or no edge) are unequal.
 Re-rewritten by Pavel Krivitsky to make compression fast. ;)
 *****************/

/* *** don't forget tail -> head, and so this function accepts
   tails before heads now */

SEXP MPLE_wrapper(SEXP stateR,
                  // MPLE settings
                  SEXP wl,
		  SEXP maxNumDyads){
  GetRNGstate(); /* Necessary for R random number generator */
  ErgmState *s = ErgmStateInit(stateR, ERGM_STATE_NO_INIT_PROP);

  Model *m = s->m;

  double *tmp = REAL(wl);
  RLEBDM1D wlm = unpack_RLEBDM1D(&tmp);

  StoreDVecMapENE *covfreq = MpleInit_hash_wl_RLE(s, &wlm, asInteger(maxNumDyads));

  // Now, unpack the hash table.
  unsigned int ntypes = kh_size(covfreq), nstats = m->n_stats;
  // R expects column-major order, so give it that.
  SEXP responsemat = PROTECT(allocMatrix(INTSXP, 2, ntypes));
  SEXP covmat = PROTECT(allocMatrix(REALSXP, nstats, ntypes));
  int *y = INTEGER(responsemat);
  double *x = REAL(covmat);

  // Save the results.
  double *k;
  ENE v;
  unsigned int i = 0;
  kh_foreach(covfreq, k, v, {
      memcpy(x + i*nstats, k, nstats*sizeof(double));
      // R's glm() expects (successes,failures)
      y[i*2u+1] = v.nonedges;
      y[i*2u] = v.edges;
      i++;
    });

  MPLE_workspace_free();

  const char *outn[] = {"y", "x", ""};
  SEXP outl = PROTECT(mkNamed(VECSXP, outn));
  SET_VECTOR_ELT(outl, 0, responsemat);
  SET_VECTOR_ELT(outl, 1, covmat);

  ErgmStateDestroy(s);
  PutRNGstate(); /* Must be called after GetRNGstate before returning to R */
  UNPROTECT(3);
  return outl;
}


static inline void insCovMatRow(StoreDVecMapENE *h, double *pred, int response){
  size_t nstat = h->l;
  int ret;

  khiter_t pos = kh_put(DVecMapENE, h, pred, &ret);
  if(ret){ // New element inserted:
    // Copy and replace the key, since it'll get overwritten later.
    double *newpred = MPLE_workspace_push(Calloc(nstat, double));
    memcpy(newpred, pred, nstat*sizeof(double));
    kh_key(h, pos) = newpred;
    kh_val(h, pos) = (ENE){0,0}; // Initialize the counters, just in case.
  }

  if(response) kh_val(h, pos).edges++;
  else kh_val(h, pos).nonedges++;
}

// Euclid's Algorithm to compute Greatest Common Divisor of a and b.
Dyad gcd(Dyad a, Dyad b){
  if(b==0) return a;
  else return gcd(b, a%b);
}

StoreDVecMapENE *MpleInit_hash_wl_RLE(ErgmState *s, RLEBDM1D *wl, Edge maxNumDyads){
  Network *nwp = s->nwp;
  Model *m = s->m;

  // Initialise the output.
  StoreDVecMapENE *covfreq = MPLE_covfreq = kh_init(DVecMapENE);
  covfreq->l = m->n_stats;

  // Number of free dyads.
  Dyad dc = wl->ndyads;
  
  // Find a number relatively prime with the free dyad count:
  Edge step = MAX(N_NODES/3,2);
  while(gcd(dc,step)!=1) step++;

  { // Start a scope for loop variables.
    Vertex t, h;
    GetRandRLEBDM1D_RS(&t,&h, wl); // Find a random starting point.
    Dyad d = TH2Dyad(nwp->nnodes, t,h);
    RLERun r=0;
    
    for(Dyad i = 0; i < MIN(maxNumDyads,dc); i++, d=NextRLEBDM1D(d, step, wl, &r)){
      R_CheckUserInterruptEvery(1024u, i);
      Dyad2TH(&t, &h, d, N_NODES);
      
      int response = IS_OUTEDGE(t,h);
      ChangeStats1(t, h, nwp, m, response);
      if(response){
	for(unsigned int l=0; l<m->n_stats; l++){
	  m->workspace[l] = -m->workspace[l];
	}
      }

      insCovMatRow(covfreq, m->workspace, response);
    }
  } // End scope for loop variables.

  return covfreq;
}
