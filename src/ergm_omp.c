/*  File src/ergm_omp.c in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2022 Statnet Commons
 */
#include <Rinternals.h>
#include "ergm_omp.h"

#ifdef _OPENMP
int ergm_omp_terms = 0;

SEXP set_ergm_omp_terms(SEXP x){
  x = PROTECT(coerceVector(x, INTSXP));
  ergm_omp_terms = INTEGER(x)[0];
  UNPROTECT(1);
  return(R_NilValue);
}

SEXP get_ergm_omp_terms(void){
  SEXP out = PROTECT(allocVector(INTSXP, 1));
  INTEGER(out)[0] = ergm_omp_terms;
  UNPROTECT(1);
  return(out);
}

#else // _OMP

SEXP set_ergm_omp_terms(SEXP x){
  error("The package was compiled without OpenMP.");
}

SEXP get_ergm_omp_terms(void){
  error("The package was compiled without OpenMP.");
}

#endif // _OMP
