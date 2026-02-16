/*  File src/allstatistics.c in package ergm, part of the Statnet suite of
 *  packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */
#include "MPLE.h"
#include "ergm_changestat.h"
#include "ergm_state.h"
#include "ergm_util.h"

KHASH_INIT(DVecMapUInt, double*, unsigned int, true, kh_DVec_hash_func, kh_DVec_hash_equal, size_t l;)
typedef khash_t(DVecMapUInt) StoreDVecMapUInt;

static kvec_t(double*) allstats_workspace = kv_blank;
static StoreDVecMapUInt *allstats_freq = NULL;

// TODO: Consider preallocating space for these vectors parcelling them out.
static inline double *allstats_workspace_push(double *ptr){
  kv_push(double*, allstats_workspace, ptr);
  return ptr;
}

SEXP allstats_workspace_free(void){
  if(allstats_freq){
    kh_destroy(DVecMapUInt, allstats_freq);
    allstats_freq = NULL;
  }

  if(kv_size(allstats_workspace))
    for(unsigned int i = 0; i < kv_size(allstats_workspace); i++)
      R_Free(kv_A(allstats_workspace, i));
  kv_destroy(allstats_workspace);

  return R_NilValue;
}

void RecurseOffOn(ErgmState *s, RLEBDM1D *wl, Dyad d, RLERun hint, double *cumulativeStats, StoreDVecMapUInt *statfreq);
        
/* *****************
 void AllStatistics

 Wrapper for a call from R.  Based on MPLE_wrapper but much simpler.  
 Produces matrix of network statistics for an arbitrary statnet model
 by an algorithm that starts with the network passed in the
 formula, then recursively toggling each edge two times so that every 
 possible network is visited.
 *****************/

SEXP AllStatistics(SEXP stateR,
                   // Allstats settings
                   SEXP wl){

  /* Step 1:  Initialize empty network and initialize model */
  GetRNGstate(); /* Necessary for R random number generator */
  ErgmState *s = ErgmStateInit(stateR, ERGM_STATE_NO_INIT_PROP);

  Model *m = s->m;
  double *tmp = REAL(wl);
  RLEBDM1D wlm = unpack_RLEBDM1D(&tmp);

  StoreDVecMapUInt *statfreq = allstats_freq = kh_init(DVecMapUInt);
  statfreq->l = m->n_stats;

  /* Step 3:  Initialize values of mtp->dstats so they point to the correct
  spots in the newRow vector.  These values will never change. */
  double *cumulativeStats = allstats_workspace_push(R_Calloc(m->n_stats,double));

  /* Step 4:  Begin recursion */
  RecurseOffOn(s, &wlm, FirstRLEBDM1D(&wlm), 1, cumulativeStats, statfreq);

  /* Step 5:  Store results and return */
  unsigned int ntypes = kh_size(statfreq), nstats = m->n_stats;
  SEXP statsR = PROTECT(allocMatrix(REALSXP, nstats, ntypes));
  SEXP weightsR = PROTECT(allocVector(INTSXP, ntypes));

  double *stats = REAL(statsR);
  int *weights = INTEGER(weightsR);

  double *k;
  unsigned int v;
  unsigned int i = 0;
  kh_foreach(statfreq, k, v, {
      memcpy(stats + i*nstats, k, nstats*sizeof(double));
      weights[i] = v;
      i++;
    });
  const char *outn[] = {"stats", "weights", ""};
  SEXP outl = PROTECT(mkNamed(VECSXP, outn));
  SET_VECTOR_ELT(outl, 0, statsR);
  SET_VECTOR_ELT(outl, 1, weightsR);


  ErgmStateDestroy(s);
  PutRNGstate(); /* Must be called after GetRNGstate before returning to R */
  UNPROTECT(3);
  return outl;
}


static inline void InsNetStatRow(StoreDVecMapUInt *h, double *newRow){
  size_t nstat = h->l;
  kh_put_code ret;
  khiter_t pos = kh_put(DVecMapUInt, h, newRow, &ret);

  if(ret != kh_put_present){ // New element inserted:
    // Copy and replace the key, since it'll get overwritten later.
    double *newstat = allstats_workspace_push(R_Calloc(nstat, double));
    memcpy(newstat, newRow, nstat*sizeof(double));
    kh_key(h, pos) = newstat;
    kh_val(h, pos) = 1;
  }else kh_val(h, pos)++;
}


static unsigned int interrupt_steps = 0; // It's OK if this overflows.

void RecurseOffOn(ErgmState *s,
                  RLEBDM1D *wl,
                  Dyad d,
                  RLERun hint,
                  double *cumulativeStats,
                  StoreDVecMapUInt *statfreq) {

  R_CheckUserInterruptEvery(1024u, interrupt_steps++);
  Network *nwp = s->nwp;
  Model *m = s->m;

  /* Current dyad */
  Vertex t, h;
  Dyad2TH(&t, &h, d, nwp->nnodes);

  /* Next dyad */
  RLERun nhint = hint;
  Dyad nd = NextRLEBDM1D(d, 1, wl, &nhint);

  /* Loop through twice for each dyad: Once for edge and once for no edge */
  for (int i=0; i<2; i++) {
    /* Recurse if currentnodes+1 is not yet nodelistlength */
    if (nd != FirstRLEBDM1D(wl)) { /* We've wrapped around. */
      RecurseOffOn(s, wl, nd, nhint, cumulativeStats, statfreq);
    } else { /* Add newRow of statistics to hash table */
      InsNetStatRow(statfreq, cumulativeStats);
    }

    /* Calculate the change statistic(s) associated with toggling the 
       dyad. */
    Rboolean edgestate = IS_OUTEDGE(t, h);
    ChangeStats1(t, h, nwp, m, edgestate);
    addonto(cumulativeStats, m->workspace, m->n_stats);

    /* Now toggle the dyad so it's ready for the next pass */
    /* Inform u_* functions that the network is about to change. */
    ToggleKnownEdge(t, h, nwp, edgestate);
  }
}
