#include "ergm_changestat_operator.h"
#include "ergm_changestat_auxnet.h"
#include "ergm_util.h"

Model *unpack_Model_as_double(double **x){
  int n_terms = *((*x)++);
  char *fnames = (char *) unpack_str_as_double(x);
  char *snames = (char *) unpack_str_as_double(x);
  Model *m = ModelInitialize(fnames, snames, x, n_terms);
  Free(fnames);
  Free(snames);
  return m;  
}
