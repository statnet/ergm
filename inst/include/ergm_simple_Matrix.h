/*  File inst/include/ergm_simple_Matrix.h in package ergm, part of the Statnet
 *  suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#ifndef _ERGM_SIMPLE_MATRIX_H_
#define _ERGM_SIMPLE_MATRIX_H_

#include "ergm_Rutil.h"

typedef struct {
  unsigned int *i, *p, dim[2];
  double *x;
} dgCMatrix_double;

static inline dgCMatrix_double *dgCMatrix_double_SEXP(SEXP from) {
  dgCMatrix_double *to = R_Calloc(1, dgCMatrix_double);
  to->i = (unsigned int *) INTEGER(R_do_slot(from, Rf_install("i")));
  to->p = (unsigned int *) INTEGER(R_do_slot(from, Rf_install("p")));
  to->x = REAL(R_do_slot(from, Rf_install("x")));
  memcpy(to->dim, INTEGER(R_do_slot(from, Rf_install("Dim"))), sizeof(int) * 2);
  return to;
}

#endif /* _ERGM_SIMPLE_MATRIX_H_ */
