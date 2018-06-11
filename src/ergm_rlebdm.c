/*  File src/ergm_rlebdm.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2018 Statnet Commons
 */
#include <R.h>
#include "ergm_rlebdm.h"

void PrintRLEBDM1D(const RLEBDM1D *m){
  Rprintf("Note: the following matrix is printed transposed:\n");
  Dyad d = 1, dmax = m->n*m->n;
  for(RLERun r = 1; r <= m->nruns; r++){
    while(d < m->starts[r-1]){
      Rprintf(".");
      if(d%m->n == 0) Rprintf("\n");
      d++;
    }
    Dyad rend =  m->starts[r-1]+m->cumlens[r]-m->cumlens[r-1];
    while(d < rend){
      Rprintf("*");
      if(d%m->n == 0) Rprintf("\n");
      d++;
    }
  }
  while(d <= dmax){
    Rprintf(".");
    if(d%m->n == 0) Rprintf("\n");
    d++;
  }
}
