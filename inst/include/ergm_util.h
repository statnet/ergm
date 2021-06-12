/*  File inst/include/ergm_util.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2021 Statnet Commons
 */
#ifndef _ERGM_UTIL_H_
#define _ERGM_UTIL_H_

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

#endif 
