/*  File inst/include/ergm_util.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2023 Statnet Commons
 */
#ifndef _ERGM_UTIL_H_
#define _ERGM_UTIL_H_

#include <R.h>

/* Calculate the dot product between two vectors. */
static inline double dotprod(double *x, double *y, unsigned int n){
  double out = 0;
  for(unsigned int i = 0; i < n; i++, x++, y++){
    out += *x * *y;
  }
  return out;
}

/* Add y to x elementwise in place. */
static inline double *addonto(double *x, double *y, unsigned int n){
  for(unsigned int i = 0; i < n; i++, x++, y++){
    *x += *y;
  }
  return x;
}

/* Pretty-print a numeric vector. */
static inline void print_vector(const char *name, double *x, unsigned int n){
  if(name) Rprintf("%s: ", name);
  Rprintf("( ");
  for(unsigned int i=0; i<n; i++)
    Rprintf("% f ", x[i]);
  Rprintf(")");
  if(name) Rprintf("\n");
}

/* Pretty-print a matrix in column-major order. */
static inline void print_matrix(const char *name, double *x, unsigned int n, unsigned int m){
  if(name) Rprintf("%s:\n", name);
  for(unsigned int i=0; i<n; i++){
    Rprintf("[ ");
    for(unsigned int j=0; j<m; j++)
      Rprintf("% f ", x[i + j*n]);
    Rprintf("]\n");
  }
}

#endif 
