#include <R.h>
#include "ergm_rlebdm.h"

void PrintBoolRLESqMatrixD(const BoolRLESqMatrixD *m){
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
