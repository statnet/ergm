#include "changestat_operator.h"

unsigned char *unpack_strasdouble(double **x){
  unsigned int l = (*x)[0];
  unsigned char *s = (unsigned char *) Calloc(l+1, char);
  for(unsigned int i=0; i<l; i++){
    s[i] = (unsigned char) (unsigned int) (*x)[i+1];
  }
  s[l] = (unsigned char) 0;
  (*x)+=l+1;
  return s;
}

Model *unpack_Modelasdouble(double **x){
  int n_terms = *((*x)++);
  char *fnames = (char *) unpack_strasdouble(x);
  char *snames = (char *) unpack_strasdouble(x);
  Model *m = ModelInitialize(fnames, snames, x, n_terms);
  Free(fnames);
  Free(snames);
  return m;  
}
