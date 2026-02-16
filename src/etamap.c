/*  File src/etamap.c in package ergm, part of the Statnet suite of packages
 *  for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */
#include "ergm_etamap.h"
#include <Rversion.h>

/* UINT_MAX's are there to segfault as soon as possible if empty from
   and to vectors are actually used rather than return misleading
   results. */
#define SETUP_CALL(fun)                                                 \
  SEXP cm = VECTOR_ELT(curved, i);                                      \
  SEXP toR = getListElement(cm, "to");                                  \
  unsigned int nto = length(toR);                                       \
  unsigned int to = nto ? INTEGER(toR)[0] : UINT_MAX;                   \
  SEXP fromR = getListElement(cm, "from");                              \
  unsigned int nfrom = length(fromR);                                   \
  unsigned int from = nfrom ? INTEGER(fromR)[0] : UINT_MAX;             \
  SEXP cov = getListElement(cm, "cov");                                 \
  SEXP fun = getListElement(cm, #fun);                                  \
                                                                        \
  SEXP pos = call, arg;                                                 \
  SETCAR(pos, fun); pos = CDR(pos);                                     \
  SETCAR(pos, (arg = allocVector(REALSXP, nfrom))); pos = CDR(pos); /* Don't need to PROTECT the vector this way. */ \
  if(nfrom) memcpy(REAL(arg), theta1+from, nfrom*sizeof(double));       \
  SETCAR(pos, ScalarInteger(nto)); pos = CDR(pos);                      \
  SETCAR(pos, cov);


/*
   A local implementation of allocLang(), recommended by Writing R
   Extensions.

   TODO: Delete in mid 2026.
*/
#if R_VERSION < R_Version(4, 4, 1)
static inline SEXP allocLang(int n)
{
    if (n > 0)
	  return LCONS(R_NilValue, allocList(n - 1));
    else
	  return R_NilValue;
}
#endif
/* End local implementation of allocLang(). */

SEXP ergm_eta_wrapper(SEXP thetaR, SEXP etamap){
  unsigned int neta = asInteger(getListElement(etamap, "etalength"));
  SEXP etaR = PROTECT(allocVector(REALSXP, neta));
  double *eta = REAL(etaR);
  double *theta = REAL(thetaR);

  ergm_eta(theta, etamap, eta);

  UNPROTECT(1);
  return etaR;
}

void ergm_eta(double *theta, SEXP etamap, double *eta){
  // For convenience, eta and theta indexed from 1:
  double *eta1 = eta-1, *theta1 = theta-1;

  SEXP ecR = getListElement(etamap, "canonical");
  unsigned int ntheta = length(ecR);
  unsigned int *ec1 = (unsigned int *) INTEGER(ecR) - 1; // Indexed from 1.
  /* Canonical parameters: */
  for(unsigned int i = 1; i <= ntheta; i++)
    if(ec1[i]) eta1[ec1[i]] = theta1[i];

  /* Curved parameters: */
  SEXP curved = getListElement(etamap, "curved");
  unsigned int ncurved = length(curved);

  if(ncurved){
    SEXP call = PROTECT(allocLang(4));

    for(unsigned int i = 0; i < ncurved; i++){
      SETUP_CALL(map);
      memcpy(eta1+to, REAL(eval(call, R_EmptyEnv)), nto*sizeof(double));
    }

    UNPROTECT(1);
  }
}

SEXP ergm_etagrad_wrapper(SEXP thetaR, SEXP etamap){
  unsigned int neta = asInteger(getListElement(etamap, "etalength"));
  SEXP etagradR = PROTECT(allocMatrix(REALSXP, length(thetaR), neta));
  double *etagrad = REAL(etagradR);
  double *theta = REAL(thetaR);

  ergm_etagrad(theta, etamap, etagrad);

  UNPROTECT(1);
  return etagradR;
}

void ergm_etagrad(double *theta, SEXP etamap, double *etagrad){
  // For convenience, eta and theta indexed from 1:
  double *theta1 = theta-1;

  SEXP ecR = getListElement(etamap, "canonical");
  unsigned int ntheta = length(ecR);
  unsigned int neta = asInteger(getListElement(etamap, "etalength"));
  memset(etagrad, 0, ntheta*neta*sizeof(double));
  double *etagrad1 = etagrad - 1 - ntheta; // Indexed from 1 (in both dimensions).

  unsigned int *ec1 = (unsigned int *) INTEGER(ecR) - 1; // Indexed from 1.
  /* Canonical parameters: */
  for(unsigned int i = 1; i <= ntheta; i++)
    if(ec1[i]) etagrad1[i+ec1[i]*ntheta] = 1;

  /* Curved parameters: */
  SEXP curved = getListElement(etamap, "curved");
  unsigned int ncurved = length(curved);

  if(ncurved){
    SEXP call = PROTECT(allocLang(4));

    for(unsigned int i = 0; i < ncurved; i++){
      SETUP_CALL(gradient);
      if(nfrom == 0) continue;
      double *g = REAL(eval(call, R_EmptyEnv));
      double *dest = etagrad1+from+to*ntheta;
      for(unsigned int j=0; j<nto; j++, dest+=ntheta, g+=nfrom)
        memcpy(dest, g, nfrom*sizeof(double));
    }

    UNPROTECT(1);
  }
}

SEXP ergm_etagradmult_wrapper(SEXP thetaR, SEXP v, SEXP etamap){
  unsigned int neta = asInteger(getListElement(etamap, "etalength"));
  unsigned int nv = isMatrix(v) ? ncols(v) : 1;
  if(neta!=(isMatrix(v) ? nrows(v) : length(v)))
    error("Non-conforming matrix multiply: grad(eta) %%*%% v.\ngrad(eta) has %u columns, but v has %u rows.",
          neta, isMatrix(v) ? nrows(v) : length(v));

  SEXP ansR = PROTECT(allocMatrix(REALSXP, length(thetaR), nv));
  double *ans = REAL(ansR);
  double *theta = REAL(thetaR);

  ergm_etagradmult(theta, REAL(v), nv, etamap, ans);

  UNPROTECT(1);
  return ansR;
}

void ergm_etagradmult(double *theta, double *v, unsigned int nv, SEXP etamap, double *ans){
  // For convenience, eta and theta indexed from 1:
  double *theta1 = theta-1;

  SEXP ecR = getListElement(etamap, "canonical");
  unsigned int ntheta = length(ecR);
  unsigned int neta = asInteger(getListElement(etamap, "etalength"));
  memset(ans, 0, ntheta*nv*sizeof(double));
  double *ans1 = ans - 1 - ntheta; // Indexed from 1 (in both dimensions).
  double *v1 = v - 1 - neta; // Indexed from 1 (in both dimensions).

  unsigned int *ec1 = (unsigned int *) INTEGER(ecR) - 1; // Indexed from 1.
  /* Canonical parameters: */
  for(unsigned int i = 1; i <= ntheta; i++)
    if(ec1[i]){
      // If these were in row-major order, we could just memcpy, but alas...
      double *dest = ans1+i+ntheta, *src = v1 + ec1[i] + neta;
      for(unsigned int j = 1; j <= nv; j++, dest+=ntheta, src+=neta)
        *dest = *src;
    }

  /* Curved parameters: */
  SEXP curved = getListElement(etamap, "curved");
  unsigned int ncurved = length(curved);

  if(ncurved){
    SEXP call = PROTECT(allocLang(4));

    for(unsigned int i = 0; i < ncurved; i++){
      SETUP_CALL(gradient);
      if(nfrom == 0) continue;
      double *g = REAL(eval(call, R_EmptyEnv));
      double *g1 = g - 1 - nfrom;
      for(unsigned int j=1; j<=nfrom; j++){
        for(unsigned int vc=1; vc<=nv; vc++){
          double dest = 0;
          for(unsigned int k=1; k<=nto; k++)
            dest += g1[j + k*nfrom] * v1[to + k-1 + vc*neta];
          ans1[from+j-1+vc*ntheta] = dest;
        }
      }
    }

    UNPROTECT(1);
  }
}
