/*  File src/ergm_Rutil.c in package ergm, part of the Statnet suite of
 *  packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#include <R.h>
#include "ergm_Rutil.h"
#include "ergm_constants.h"

int GetBuiltABIVersion_ergm(void){
  return ERGM_ABI_VERSION;
}

SEXP GetBuiltABIVersion_wrapper(SEXP client, SEXP lib){
  const char *sn = FIRSTCHAR(client), *ln = FIRSTCHAR(lib);
  int lfn = strlen("GetBuiltABIVersion_") + strlen(ln);
  char *fn = R_Calloc(lfn + 1, char);

  snprintf(fn, lfn + 1, "GetBuiltABIVersion_%s", ln);
  for(unsigned int i = 0; i < lfn; i++) if(fn[i] == '.') fn[i] = '_';

  int (*f)(void) = (int (*)(void)) R_FindSymbol(fn, sn, NULL);

  R_Free(fn);

  return f ?  ScalarInteger(f()) : R_NilValue;
}


SEXP mat_by_coef(SEXP A, SEXP x) {
  // Ensure A is a numeric matrix
  if (!isReal(A) || !isMatrix(A)) {
    error("A must be a numeric matrix.");
  }
  // Ensure x is a numeric vector
  if (!isReal(x) || !isVector(x)) {
    error("x must be a numeric vector.");
  }

  SEXP dimA = getAttrib(A, R_DimSymbol);
  int m = INTEGER(dimA)[0];
  int n = INTEGER(dimA)[1];

  if (LENGTH(x) != n) {
    error("Length of x (%d) does not match number of columns in A (%d).", LENGTH(x), n);
  }

  const double *a = REAL(A);
  const double *xx = REAL(x);

  SEXP ans = PROTECT(allocVector(REALSXP, m));
  double *rans = REAL(ans);

  for (int i = 0; i < m; i++) {
    double sum = 0.0;
    for (int j = 0; j < n; j++) {
      double aij = a[i + j * m];
      double xj  = xx[j];

      double prod;
      if ((aij == 0.0 && (isinf(xj) || isnan(xj))) ||
          (xj == 0.0 && (isinf(aij) || isnan(aij)))) {
        prod = 0.0;
      } else {
        prod = aij * xj;
      }
      sum += prod;
    }
    rans[i] = sum;
  }

  UNPROTECT(1);
  return ans;
}
