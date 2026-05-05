#include <R.h>
#include <Rinternals.h>

SEXP ar_ols_stats(SEXP Xl, SEXP lagSEXP) {
  const int L = asInteger(lagSEXP);
  const int nlist = LENGTH(Xl);

  SEXP X0 = VECTOR_ELT(Xl, 0);
  const int p = INTEGER(getAttrib(X0, R_DimSymbol))[1];
  const int K = p * L;

  /* Lagged cross-products */
  SEXP XtX = PROTECT(allocMatrix(REALSXP, K, K));   /* 1 */
  SEXP XtY = PROTECT(allocMatrix(REALSXP, K, p));   /* 2 */

  /* Conditional-on-L response moments */
  SEXP YtY_L = PROTECT(allocMatrix(REALSXP, p, p)); /* 3 */
  SEXP oneY_L = PROTECT(allocMatrix(REALSXP, 1, p));/* 4 */
  SEXP one1_L = PROTECT(ScalarInteger(0));          /* 5 */

  double *XtXr = REAL(XtX);
  double *XtYr = REAL(XtY);
  double *YtYr = REAL(YtY_L);
  double *oneYr = REAL(oneY_L);

  memset(XtXr,  0, sizeof(double) * (size_t)K * (size_t)K);
  memset(XtYr,  0, sizeof(double) * (size_t)K * (size_t)p);
  memset(YtYr,  0, sizeof(double) * (size_t)p * (size_t)p);
  memset(oneYr, 0, sizeof(double) * (size_t)p);

  int nobs_L = 0;

  for (int c = 0; c < nlist; c++) {
    SEXP X = VECTOR_ELT(Xl, c);
    double *Xr = REAL(X);
    const int T = INTEGER(getAttrib(X, R_DimSymbol))[0];

    /* Conditional response moments: t = L+1 ... T */
    if (T > L) {
      const int n = T - L;
      nobs_L += n;

      for (int a = 0; a < p; a++) {
        double *xi = Xr + T*a + L, acc = 0;
        for (int t = 0; t < n; t++) acc += *xi++;
        oneYr[a] += acc;

        for (int b = 0; b < p; b++) {
          double *xi = Xr + T*a + L, *xj = Xr + T*b + L, acc = 0;
          for (int t = 0; t < n; t++) acc += (*xi++) * (*xj++);
          YtYr[a + p*b] += acc;
        }
      }
    }

    const int Lmax = L >= T ? T-1 : L;
    for (int i = 1; i <= Lmax; i++) {
      const int n = T - i;
      /* XtY */
      for (int a = 0; a < p; a++)
        for (int b = 0; b < p; b++) {
          double *xi = Xr + T*a, *y = Xr + T*b + i, acc = 0;
          for (int t = 0; t < n; t++) acc += (*xi++) * (*y++);

          XtYr[(i-1)*p + a + K*b] += acc;
        }

      /* XtX */
      for (int j = 1; j <= i; j++)
        for (int a = 0; a < p; a++)
          for (int b = 0; b < p; b++) {
            double *xi = Xr + T*a, *xj = Xr + T*b + (i - j), acc = 0;
            for (int t = 0; t < n; t++) acc += (*xi++) * (*xj++);

            /* lower triangle */
            XtXr[(i-1)*p + a + K*((j-1)*p + b)] += acc;

            /* mirrored upper triangle */
            if (i != j) XtXr[(j-1)*p + b + K*((i-1)*p + a)] += acc;
          }
    }
  }

  INTEGER(one1_L)[0] = nobs_L;

  static const char *names[] = {
    "XtX", "XtY", "YtY_L", "oneY_L", "one1_L", ""
  };

  SEXP out = PROTECT(Rf_mkNamed(VECSXP, names)); /* 6 */

  SET_VECTOR_ELT(out, 0, XtX);
  SET_VECTOR_ELT(out, 1, XtY);
  SET_VECTOR_ELT(out, 2, YtY_L);
  SET_VECTOR_ELT(out, 3, oneY_L);
  SET_VECTOR_ELT(out, 4, one1_L);

  UNPROTECT(6);
  return out;
}
