/*  File src/wtnetstats.c in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2022 Statnet Commons
 */
#include "ergm_wtstate.h"
/*****************
 void network_stats_wrapper

 Wrapper for a call from R.  Return the statistics when of the observed graph.
*****************/

SEXP wt_network_stats_wrapper(SEXP stateR){
  GetRNGstate();  /* R function enabling uniform RNG */
  WtErgmState *s = WtErgmStateInit(stateR,
                                   ERGM_STATE_EMPTY_NET | ERGM_STATE_NO_INIT_S | ERGM_STATE_NO_INIT_PROP);

  WtModel *m = s->m;

  SEXP stats = PROTECT(allocVector(REALSXP, m->n_stats));
  m->workspace = REAL(stats);

  SEXP elR = getListElement(stateR, "el");
  WtSummStats(length(VECTOR_ELT(elR, 0)),
              (Vertex*) INTEGER(VECTOR_ELT(elR, 0)),
              (Vertex*) INTEGER(VECTOR_ELT(elR, 1)),
              REAL(VECTOR_ELT(elR, 2)),
              s->nwp, m);

  WtErgmStateDestroy(s);
  PutRNGstate();
  UNPROTECT(1);
  return stats;
}
