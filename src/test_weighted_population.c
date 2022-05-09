/*  File src/test_weighted_population.c in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2022 Statnet Commons
 */

#include "ergm_weighted_population.h"
#include <Rinternals.h>

SEXP test_weighted_population(SEXP weights, SEXP ndraws, SEXP type) {
  GetRNGstate();  /* R function enabling uniform RNG */
  
  WtPop *wtp = WtPopInitialize(length(weights), REAL(weights), CHAR(asChar(type))[0]);
  
  int n = asInteger(ndraws);
  
  SEXP draws = PROTECT(allocVector(INTSXP, n));
  int *d = INTEGER(draws);
  memset(d, 0, n*sizeof(int));
  
  for(int i = 0; i < n; i++) {
    d[i] = WtPopGetRand(wtp);
  }
  
  WtPopDestroy(wtp);

  PutRNGstate();  /* Disable RNG before returning */
  
  UNPROTECT(1);
  return draws;
}
