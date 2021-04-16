#ifndef _ERGM_HASH_EDGELIST_H_
#define _ERGM_HASH_EDGELIST_H_

#include "ergm_unsorted_edgelist.h"
#include "ergm_dyad_hashmap.h"

// provides O(1) sampling, hash table time searching,
// and hash table time + O(1) inserting and deleting

typedef struct {
  UnsrtEL *list;
  StoreDyadMapUInt *hash;
} HashEL;

static inline HashEL *HashELInitialize(unsigned int nedges, Vertex *tails, Vertex *heads, Rboolean copy, Rboolean directed) {
  HashEL *hash = Calloc(1, HashEL);

  hash->list = UnsrtELInitialize(nedges, tails, heads, copy);

  hash->hash = kh_init(DyadMapUInt);
  hash->hash->directed = directed;
  
  if(nedges > 0) {
    kh_resize(DyadMapUInt, hash->hash, 2*(nedges + 1));
    
    for(unsigned int i = 0; i < nedges; i++) {
      kh_set(DyadMapUInt, hash->hash, THKey(hash->hash, tails[i], heads[i]), i + 1);
    }
  }
  
  return hash;
}

static inline void HashELDestroy(HashEL *hash) {
  kh_destroy(DyadMapUInt, hash->hash);
  UnsrtELDestroy(hash->list);
  Free(hash);
}

static inline void HashELClear(HashEL *hash) {
  UnsrtELClear(hash->list);
  kh_clear(DyadMapUInt, hash->hash);
}

static inline void HashELGetRand(Vertex *tail, Vertex *head, HashEL *hash) {
  UnsrtELGetRand(tail, head, hash->list);
}

static inline void HashELInsert(Vertex tail, Vertex head, HashEL *hash) {
  UnsrtELInsert(tail, head, hash->list);
  kh_set(DyadMapUInt, hash->hash, THKey(hash->hash, tail, head), hash->list->nedges);
}

static inline void HashELDelete(Vertex tail, Vertex head, HashEL *hash) {   
  khint_t i = kh_get(DyadMapUInt, hash->hash, THKey(hash->hash, tail, head));    
  unsigned int index = kh_value(hash->hash, i);
  kh_del(DyadMapUInt, hash->hash, i);

  if(index < hash->list->nedges) {
    kh_set(DyadMapUInt, hash->hash, THKey(hash->hash, hash->list->tails[hash->list->nedges], hash->list->heads[hash->list->nedges]), index);
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
  return kh_getval(DyadMapUInt, hash->hash, THKey(hash->hash, tail, head), 0);
}

#endif // _ERGM_HASH_EDGELIST_H_
