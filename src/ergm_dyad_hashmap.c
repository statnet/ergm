#include "ergm_dyad_hashmap.h"

/* Print the contents of a khash table  mapping dyads to unsigned
   integers. Useful for debugging. */
void PrintDyadMapUInt(StoreDyadMapUInt *h){
  for(khiter_t i = kh_begin(h); i!=kh_end(h); ++i){
    if(kh_exist(h, i)){
      struct TailHead k = kh_key(h, i);
      unsigned int v = kh_val(h, i);
      Rprintf("(%d,%d)->%u\n",k.tail,k.head,v);
    }
  }
}
