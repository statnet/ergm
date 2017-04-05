#include "changestat_operator.h"
#include "ergm_util.h"

Model *unpack_Modelasdouble(double **x){
  int n_terms = *((*x)++);
  char *fnames = (char *) unpack_strasdouble(x);
  char *snames = (char *) unpack_strasdouble(x);
  Model *m = ModelInitialize(fnames, snames, x, n_terms);
  Free(fnames);
  Free(snames);
  return m;  
}
