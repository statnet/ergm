/*  File inst/include/ergm_dyad_hashmap.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2021 Statnet Commons
 */
#ifndef _ERGM_DYAD_HASHMAP_H_
#define _ERGM_DYAD_HASHMAP_H_

#include <R.h>
#include "ergm_edgetree_types.h"

/* Specify allocators. */
#define kcalloc(N,Z) R_chk_calloc(N,Z)
#define kmalloc(Z) R_chk_calloc(Z,1)
#define krealloc(P,Z) R_chk_realloc(P,Z)
#define kfree(P) R_chk_free(P)
#include "ergm_khash.h"

/* Data structure to represent a dyad that can serve as a key to the hash. */
typedef struct TailHead_s{
  Vertex tail, head;
} TailHead;

/* Helper macros to construct TailHeads with the correct ordering. */
#define TH(T,H) ((TailHead){.tail=(T),.head=(H)})

/* Hash and comparison functions designed for tail-head pairs. */
// The of macro-ing this due to Bob Jenkins.
#define ROT_INT(x,k) (((x)<<(k)) | ((x)>>(32-(k))))
// This combination of rotations has been taken from
// http://faculty.otterbein.edu/psanderson/csc205/notes/lecture17.html
// and is used in Java 5ish. It has been tested extensively.
static inline unsigned int kh_scramble_int(unsigned int a){
  a += ~(a << 9);
  a ^= (a >> 14);
  a += ~(a << 4);
  a ^= (a >> 10);
  return a;
}
/* #define kh_vertexvertex_hash_func(key) (khint32_t)(kh_scramble_int(ROT_INT((key).tail,16) ^ (key).head)) */
#define kh_vertexvertex_hash_func(key) (khint32_t) (((key).tail<(key).head || h->directed) ? ((key).tail + (key).head*0xd7d4eb2du) : ((key).head + (key).tail*0xd7d4eb2du))
#define kh_vertexvertex_hash_equal(a,b) ((a.tail==b.tail && a.head==b.head) || (!h->directed && a.tail==b.head && a.head==b.tail))

/* Predefined khash type for mapping dyads onto unsigned ints. */
KHASH_INIT(DyadMapUInt, TailHead, unsigned int, TRUE, kh_vertexvertex_hash_func, kh_vertexvertex_hash_equal, Rboolean directed;)
typedef khash_t(DyadMapUInt) StoreDyadMapUInt;

/* Accessors, modifiers, and incrementors. */
#define GETDMUI(tail, head, hashmap)(kh_getval(DyadMapUInt, hashmap, TH(tail,head), 0))
#define SETDMUI(tail, head, v, hashmap) {if(v==0) kh_unset(DyadMapUInt, hashmap, TH(tail,head)); else kh_set(DyadMapUInt, hashmap, TH(tail,head), v)}
#define HASDMUI(tail, head, hashmap)(kh_get(DyadMapUInt, hashmap, TH(tail,head))!=kh_none)
#define DELDMUI(tail, head, hashmap)(kh_unset(DyadMapUInt, hashmap, TH(tail,head)))
#define SETDMUI0(tail, head, v, hashmap) {kh_set(DyadMapUInt, hashmap, TH(tail,head), v)}

static inline void IncDyadMapUInt(Vertex tail, Vertex head, int inc, StoreDyadMapUInt *h){
  TailHead th = TH(tail, head);
  if(inc!=0){
    khiter_t pos = kh_get(DyadMapUInt, h, th);
    unsigned int val = pos==kh_none ? 0 : kh_value(h, pos);
    val += inc;
    if(val==0){
      kh_del(DyadMapUInt, h, pos);
    }else{
      if(pos==kh_none)
        pos = kh_put(DyadMapUInt, h, th, NULL);
      kh_val(h, pos) = val;
    }
  }
}

/* Predefined khash type for mapping dyads onto signed ints. */
KHASH_INIT(DyadMapInt, TailHead, int, TRUE, kh_vertexvertex_hash_func, kh_vertexvertex_hash_equal, Rboolean directed;)
typedef khash_t(DyadMapInt) StoreDyadMapInt;

/* Accessors, modifiers, and incrementors. */
#define GETDMI(tail, head, hashmap)(kh_getval(DyadMapInt, hashmap, TH(tail,head), 0))
#define SETDMI(tail, head, v, hashmap) {if(v==0) kh_unset(DyadMapInt, hashmap, TH(tail,head)); else kh_set(DyadMapInt, hashmap, TH(tail,head), v)}
#define HASDMI(tail, head, hashmap)(kh_get(DyadMapInt, hashmap, TH(tail,head))!=kh_none)
#define DELDMI(tail, head, hashmap)(kh_unset(DyadMapInt, hashmap, TH(tail,head)))
#define SETDMI0(tail, head, v, hashmap) {kh_set(DyadMapInt, hashmap, TH(tail,head), v)}

/* Predefined khash type for dyad sets. This may or may not be faster than edgetree. */
KHASH_INIT(DyadSet, TailHead, char, FALSE, kh_vertexvertex_hash_func, kh_vertexvertex_hash_equal, Rboolean directed;)
typedef khash_t(DyadSet) StoreDyadSet;

// Toggle an element of a DyadSet.
static inline Rboolean DyadSetToggle(Vertex tail, Vertex head, StoreDyadSet *h){
  TailHead th = TH(tail, head);
  int ret;
  // Attempt insertion
  khiter_t i = kh_put(DyadSet, h, th, &ret);
  if(ret==0){
    // Already present: delete
    kh_del(DyadSet, h, i);
    return FALSE;
  }else{
    // Inserted by kh_put above
    return TRUE;
  }
}

#endif // _ERGM_DYAD_HASHMAP_H_
