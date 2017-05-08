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

SEXP get_ergm_omp_terms(){
  SEXP out = PROTECT(allocVector(INTSXP, 1));
  INTEGER(out)[0] = ergm_omp_terms;
  UNPROTECT(1);
  return(out);
}

#else // _OMP

SEXP set_ergm_omp_terms(SEXP x){
  error("The package was compiled without OpenMP.");
}

SEXP get_ergm_omp_terms(){
  error("The package was compiled without OpenMP.");
}

#endif // _OMP
