/*  File src/netstats.c.template.do_not_include_directly.h in package ergm,
 *  part of the Statnet suite of packages for network analysis,
 *  https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */

/*****************
 void NetStats_wrapper

 Wrapper for a call from R.  Return the statistics when of the observed graph.
*****************/

SEXP ETYPE(SummStats_wrapper)(SEXP stateR){
  GetRNGstate();  /* R function enabling uniform RNG */
  ETYPE(ErgmState) *s = ETYPE(ErgmStateInit)(stateR,
                                   ERGM_STATE_EMPTY_NET | ERGM_STATE_NO_INIT_S | ERGM_STATE_NO_INIT_PROP);

  ETYPE(Model) *m = s->m;

  SEXP stats = PROTECT(allocVector(REALSXP, m->n_stats));
  m->workspace = REAL(stats);

  SEXP elR = getListElement(stateR, "el");
  ETYPE(SummStats)(length(VECTOR_ELT(elR, 0)),
              (Vertex*) INTEGER(VECTOR_ELT(elR, 0)),
              (Vertex*) INTEGER(VECTOR_ELT(elR, 1)),
              IFEWT(EWTRTYPE(VECTOR_ELT(elR, 2)),)
              s->nwp, m);

  ETYPE(ErgmStateDestroy)(s);
  PutRNGstate();
  UNPROTECT(1);
  return stats;
}
