#include <Rinternals.h>
#include "ergm_omp.h"

unsigned int ergm_omp_terms = 0;

SEXP set_ergm_omp_terms(SEXP x){
  x = PROTECT(coerceVector(x, LGLSXP));
  ergm_omp_terms = LOGICAL(x)[0];
  UNPROTECT(1);
  return(R_NilValue);
}

SEXP get_ergm_omp_terms(){
  SEXP out = PROTECT(allocVector(LGLSXP, 1));
  LOGICAL(out)[0] = ergm_omp_terms;
  UNPROTECT(1);
  return(out);
}
