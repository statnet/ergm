#include "wtchangestat_operator.h"
#include "ergm_util.h"

WtModel *unpack_WtModelasdouble(double **x){
  int n_terms = *((*x)++);
  char *fnames = (char *) unpack_strasdouble(x);
  char *snames = (char *) unpack_strasdouble(x);
  WtModel *m = WtModelInitialize(fnames, snames, x, n_terms);
  Free(fnames);
  Free(snames);
  return m;  
}
