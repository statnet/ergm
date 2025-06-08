/*  File src/ergm_dyad_hashmap_utils.c in package ergm, part of the Statnet
 *  suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#include "ergm_dyad_hashmap.h"
#include "ergm_changestat.h"
#include "ergm_dyad_hashmap_utils.h"

/* Print the contents of a khash table  mapping dyads to unsigned
   integers. Useful for debugging. */
void PrintDyadMapUInt(StoreDyadMapUInt *h){
  for(khiter_t i = kh_begin(h); i!=kh_end(h); ++i){
    if(kh_exist(h, i)){
      TailHead k = kh_key(h, i);
      unsigned int v = kh_val(h, i);
      Rprintf("(%d,%d)->%u\n",k.tail,k.head,v);
    }
  }
}

/* Print the contents of a khash table  mapping dyads to unsigned
   integers. Useful for debugging. */
void PrintStrictDyadMapUInt(StoreStrictDyadMapUInt *h){
  for(khiter_t i = kh_begin(h); i!=kh_end(h); ++i){
    if(kh_exist(h, i)){
      TailHead k = kh_key(h, i);
      unsigned int v = kh_val(h, i);
      Rprintf("(%d,%d)->%u\n",k.tail,k.head,v);
    }
  }
}

/* Print the contents of a khash set of dyads. Useful for debugging. */
void PrintDyadSet(StoreDyadSet *h){
  for(khiter_t i = kh_begin(h); i!=kh_end(h); ++i){
    if(kh_exist(h, i)){
      TailHead k = kh_key(h, i);
      Rprintf("(%d,%d) ",k.tail,k.head);
    }
  }
  Rprintf("\n");
}

/* Print the contents of a khash set of dyads. Useful for debugging. */
void PrintStrictDyadSet(StoreStrictDyadSet *h){
  for(khiter_t i = kh_begin(h); i!=kh_end(h); ++i){
    if(kh_exist(h, i)){
      TailHead k = kh_key(h, i);
      Rprintf("(%d,%d) ",k.tail,k.head);
    }
  }
  Rprintf("\n");
}

/* Copy network to a khash set of dyads. */
StoreDyadSet *NetworkToDyadSet(Network *nwp){
  StoreDyadSet *h = kh_init(DyadSet);
  h->directed = DIRECTED;

  EXEC_THROUGH_NET_EDGES(tail, head, e, {
      kh_put(DyadSet, h, TH(tail,head), NULL);
    });
  return h;
}

StoreStrictDyadSet *NetworkToStrictDyadSet(Network *nwp){
  StoreStrictDyadSet *h = kh_init(StrictDyadSet);

  EXEC_THROUGH_NET_EDGES(tail, head, e, {
      kh_put(StrictDyadSet, h, TH(tail,head), NULL);
    });
  return h;
}
