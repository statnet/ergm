/*  File src/ergm_util.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2017 Statnet Commons
 */
#ifndef _ERGM_UTIL_H_
#define _ERGM_UTIL_H_
unsigned char *unpack_str_as_double(double **x);

static inline double dotprod(double *x, double *y, unsigned int n){
  double out = 0;
  for(unsigned int i = 0; i < n; i++, x++, y++){
    out += *x * *y;
  }
  return out;
}

#define TOINTSXP(x) x = PROTECT(coerceVector(x, INTSXP))
#define TOREALSXP(x) x = PROTECT(coerceVector(x, REALSXP))
#define FIRSTCHAR(x) CHAR(STRING_ELT(x, 0))

#define NWSTATE_NAMES "newnwtails", "newnwheads"
#define NWSTATE_SAVE_INTO_RLIST(nwp, outl, pos)                         \
  {                                                                     \
    SEXP newnetworktails = PROTECT(allocVector(INTSXP, EDGECOUNT(nwp)+1)); \
    SEXP newnetworkheads = PROTECT(allocVector(INTSXP, EDGECOUNT(nwp)+1)); \
    INTEGER(newnetworktails)[0]=INTEGER(newnetworkheads)[0]=            \
      EdgeTree2EdgeList((Vertex*)INTEGER(newnetworktails)+1,            \
                        (Vertex*)INTEGER(newnetworkheads)+1,            \
                        (nwp),(nwp)->nedges);                           \
    SET_VECTOR_ELT((outl), (pos), newnetworktails);                     \
    SET_VECTOR_ELT((outl), (pos)+1, newnetworkheads);                   \
    UNPROTECT(2);                                                       \
  }

#define WTNWSTATE_NAMES "newnwtails", "newnwheads", "newnwweights"
#define WTNWSTATE_SAVE_INTO_RLIST(nwp, outl, pos)                       \
  {                                                                     \
  SEXP newnetworktails = PROTECT(allocVector(INTSXP, EDGECOUNT(nwp)+1)); \
  SEXP newnetworkheads = PROTECT(allocVector(INTSXP, EDGECOUNT(nwp)+1)); \
  SEXP newnetworkweights = PROTECT(allocVector(REALSXP, EDGECOUNT(nwp)+1)); \
  INTEGER(newnetworktails)[0]=INTEGER(newnetworkheads)[0]=REAL(newnetworkweights)[0]= \
    WtEdgeTree2EdgeList((Vertex*)INTEGER(newnetworktails)+1,            \
                        (Vertex*)INTEGER(newnetworkheads)+1,            \
                        REAL(newnetworkweights)+1,                      \
                        (nwp),(nwp)->nedges);                           \
  SET_VECTOR_ELT((outl), (pos), newnetworktails);                       \
  SET_VECTOR_ELT((outl), (pos)+1, newnetworkheads);                     \
  SET_VECTOR_ELT((outl), (pos)+2, newnetworkweights);                   \
  UNPROTECT(3);                                                         \
  }

#endif 



