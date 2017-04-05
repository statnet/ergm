#include <R.h>
#include "ergm_util.h"

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

