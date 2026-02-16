/*  File inst/include/ergm_hash_edgelist.h in package ergm, part of the Statnet
 *  suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */
#ifndef _ERGM_HASH_EDGELIST_H_
#define _ERGM_HASH_EDGELIST_H_

#include "ergm_unsorted_edgelist.h"
#include "ergm_dyad_hashmap.h"

/*
   HashEL provides O(1) sampling, hash table time searching,
   and hash table time + O(1) inserting and deleting;
   unlike UnsrtEL, the HashEL data structure remains
   efficient when deleting an edge that was not just sampled
*/

typedef struct {
  UnsrtEL *list;
  StoreStrictDyadMapUInt *hash;
} HashEL;

#define HashELSize(hash) (UnsrtELSize((hash)->list))

/* Embed an existing UnsrtEL into a HashEL. */
static inline HashEL *UnsrtELIntoHashEL(UnsrtEL *el) {
  HashEL *hash = R_Calloc(1, HashEL);

  hash->list = el;

  hash->hash = kh_init(StrictDyadMapUInt);

  if(UnsrtELSize(el) > 0) {
    kh_resize(StrictDyadMapUInt, hash->hash, 2*(UnsrtELSize(el) + 1));

    for(unsigned int i = 1; i <= UnsrtELSize(el); i++) {
      kh_set(StrictDyadMapUInt, hash->hash, TH(el->tails[i], el->heads[i]), i);
    }
  }

  return hash;
}

static inline HashEL *HashELInitialize(unsigned int nedges, Vertex *tails, Vertex *heads, Rboolean copy) {
  UnsrtEL *el = UnsrtELInitialize(nedges, tails, heads, copy);
  return UnsrtELIntoHashEL(el);
}

static inline void HashELDestroy(HashEL *hash) {
  kh_destroy(StrictDyadMapUInt, hash->hash);
  UnsrtELDestroy(hash->list);
  R_Free(hash);
}

static inline void HashELClear(HashEL *hash) {
  UnsrtELClear(hash->list);
  kh_clear(StrictDyadMapUInt, hash->hash);
}

static inline void HashELGetRand(Vertex *tail, Vertex *head, HashEL *hash) {
  UnsrtELGetRand(tail, head, hash->list);
}

static inline void HashELInsert(Vertex tail, Vertex head, HashEL *hash) {
  kh_put_code r;
  khiter_t pos = kh_put(StrictDyadMapUInt, hash->hash, TH(tail, head), &r);
  if(r == kh_put_present) return; // Already in the list.

  UnsrtELInsert(tail, head, hash->list);
  kh_val(hash->hash, pos) = HashELSize(hash);
}

static inline void HashELDelete(Vertex tail, Vertex head, HashEL *hash) {
  khint_t i = kh_get(StrictDyadMapUInt, hash->hash, TH(tail, head));
  unsigned int index = kh_value(hash->hash, i);
  kh_del(StrictDyadMapUInt, hash->hash, i);

  if(index < HashELSize(hash)) {
    kh_set(StrictDyadMapUInt,
           hash->hash,
           TH(hash->list->tails[HashELSize(hash)],
              hash->list->heads[HashELSize(hash)]),
           index);
  }

  UnsrtELDeleteAt(index, hash->list);
}

static inline void HashELToggleKnown(Vertex tail, Vertex head, HashEL *hash, int edgestate) {
  if(edgestate) {
    HashELDelete(tail, head, hash);
  } else {
    HashELInsert(tail, head, hash);
  }
}

static inline unsigned int HashELSearch(Vertex tail, Vertex head, HashEL *hash) {
  return kh_getval(StrictDyadMapUInt, hash->hash, TH(tail, head), 0);
}

#endif // _ERGM_HASH_EDGELIST_H_
