/*  File src/MHproposal.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2013 Statnet Commons
 */
#ifndef _ERGM_UTIL_H_
#define _ERGM_UTIL_H_
unsigned char *unpack_str_as_double(double **x);

static inline double dotprod(double *x, double *y, unsigned int n){
  double out = 0;
  for(unsigned int i = 0; i < n; i++, x++, y++){
    out += *x * *y;
  }
  return out;
}

#endif 



