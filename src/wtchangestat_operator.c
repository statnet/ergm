#include "wtchangestat_operator.h"
#include "ergm_util.h"

WtModel *unpack_WtModel_as_double(double **x){
  int n_terms = *((*x)++);
  char *fnames = (char *) unpack_str_as_double(x);
  char *snames = (char *) unpack_str_as_double(x);
  WtModel *m = WtModelInitialize(fnames, snames, x, n_terms);
  Free(fnames);
  Free(snames);
  return m;  
}
