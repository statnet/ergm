/*  File inst/include/ergm_dyad_hashmap.h in package ergm, part of the Statnet
 *  suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#ifndef _ERGM_DYAD_HASHMAP_H_
#define _ERGM_DYAD_HASHMAP_H_

#include <R.h>
#include "ergm_edgetree_types.h"

#include "ergm_khash.h"

/* Data structure to represent a dyad that can serve as a key to the hash. */
typedef struct TailHead_s{
  Vertex tail, head;
} TailHead;

/* Helper macros to construct TailHeads with the correct ordering. */
#define TH(T,H) ((TailHead){.tail=(T),.head=(H)})
#define UTH(tail, head) tail < head ? TH(tail, head) : TH(head, tail)

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

#define kh_vertexvertex_strict_hash_func(key) (khint32_t) ((key).tail + (key).head*0xd7d4eb2du)
#define kh_vertexvertex_strict_hash_equal(a,b) (a.tail==b.tail && a.head==b.head)

/* Accessors, modifiers, and incrementors. */
#define AccessorTemplate(fname, type, valtype, th_impl)                 \
  static inline valtype GET ## fname (Vertex tail, Vertex head, Store ## type *hashmap){ \
    return kh_getval(type, hashmap, th_impl, 0);                        \
  }                                                                     \
  static inline void SET ## fname (Vertex tail, Vertex head, valtype v, Store ## type *hashmap){ \
    if(v == (valtype)0){                                                \
      kh_unset(type, hashmap, th_impl);                                 \
    }else{                                                              \
      kh_set(type, hashmap, th_impl, v);                                \
    }                                                                   \
  }                                                                     \
  static inline valtype HAS ## fname (Vertex tail, Vertex head, Store ## type *hashmap){ \
    return kh_get(type, hashmap, th_impl) != kh_none;                   \
  }                                                                     \
  static inline void DEL ## fname (Vertex tail, Vertex head, Store ## type *hashmap){ \
    kh_unset(type, hashmap, th_impl);                            \
  }                                                                     \
  static inline void SET ## fname ## 0 (Vertex tail, Vertex head, valtype v, Store ## type *hashmap){ \
    kh_set(type, hashmap, th_impl, v);                                  \
  }

#define IncDyadMapTemplate(fname, type, valtype, th_impl)               \
  static inline void fname(Vertex tail, Vertex head, int inc, Store ## type *h){ \
    TailHead th = th_impl;                                              \
    if(inc!=0){                                                         \
      khiter_t pos = kh_get(type, h, th);                               \
      valtype val = pos==kh_none ? 0 : kh_value(h, pos);                \
      val += inc;                                                       \
      if(val==0){                                                       \
        kh_del(type, h, pos);                                           \
      }else{                                                            \
        if(pos==kh_none)                                                \
          pos = kh_put(type, h, th, NULL);                              \
        kh_val(h, pos) = val;                                           \
      }                                                                 \
    }                                                                   \
  }


/* Predefined khash type for mapping dyads onto unsigned ints. */
KHASH_INIT(DyadMapUInt, TailHead, unsigned int, TRUE, kh_vertexvertex_hash_func, kh_vertexvertex_hash_equal, Rboolean directed;)
typedef khash_t(DyadMapUInt) StoreDyadMapUInt;
AccessorTemplate(DMUI, DyadMapUInt, unsigned int, TH(tail,head))
IncDyadMapTemplate(IncDyadMapUInt, DyadMapUInt, unsigned int, TH(tail, head))

KHASH_INIT(StrictDyadMapUInt, TailHead, unsigned int, TRUE, kh_vertexvertex_strict_hash_func, kh_vertexvertex_strict_hash_equal,)
typedef khash_t(StrictDyadMapUInt) StoreStrictDyadMapUInt;
AccessorTemplate(DDMUI, StrictDyadMapUInt, unsigned int, TH(tail,head))
AccessorTemplate(UDMUI, StrictDyadMapUInt, unsigned int, UTH(tail, head))
IncDyadMapTemplate(IncDDyadMapUInt, StrictDyadMapUInt, unsigned int, TH(tail, head))
IncDyadMapTemplate(IncUDyadMapUInt, StrictDyadMapUInt, unsigned int, UTH(tail, head))

/* Predefined khash type for mapping dyads onto signed ints. */
KHASH_INIT(DyadMapInt, TailHead, int, TRUE, kh_vertexvertex_hash_func, kh_vertexvertex_hash_equal, Rboolean directed;)
typedef khash_t(DyadMapInt) StoreDyadMapInt;
AccessorTemplate(DMI, DyadMapInt, int, TH(tail,head))
IncDyadMapTemplate(IncDyadMapInt, DyadMapInt, int, TH(tail, head))

KHASH_INIT(StrictDyadMapInt, TailHead, int, TRUE, kh_vertexvertex_strict_hash_func, kh_vertexvertex_strict_hash_equal,)
typedef khash_t(StrictDyadMapInt) StoreStrictDyadMapInt;
AccessorTemplate(DDMI, StrictDyadMapInt, int, TH(tail,head))
AccessorTemplate(UDMI, StrictDyadMapInt, int, UTH(tail, head))
IncDyadMapTemplate(IncDDyadMapInt, StrictDyadMapInt, int, TH(tail, head))
IncDyadMapTemplate(IncUDyadMapInt, StrictDyadMapInt, int, UTH(tail, head))

/* Predefined khash type for dyad sets. This may or may not be faster than edgetree. */

// Toggle an element of a DyadSet.
#define DyadSetToggleTemplate(fname, type, th_impl)                     \
  static inline Rboolean fname(Vertex tail, Vertex head, Store ## type *h){ \
    TailHead th = th_impl;                                              \
    kh_put_code ret;                                                    \
    khiter_t i = kh_put(type, h, th, &ret);                             \
    if(ret==0){                                                         \
      kh_del(type, h, i);                                               \
      return FALSE;                                                     \
    }else{                                                              \
      return TRUE;                                                      \
    }                                                                   \
  }

KHASH_INIT(DyadSet, TailHead, char, FALSE, kh_vertexvertex_hash_func, kh_vertexvertex_hash_equal, Rboolean directed;)
typedef khash_t(DyadSet) StoreDyadSet;
DyadSetToggleTemplate(DyadSetToggle, DyadSet, TH(tail, head))

KHASH_INIT(StrictDyadSet, TailHead, char, FALSE, kh_vertexvertex_strict_hash_func, kh_vertexvertex_strict_hash_equal,)
typedef khash_t(StrictDyadSet) StoreStrictDyadSet;
DyadSetToggleTemplate(DDyadSetToggle, StrictDyadSet, TH(tail, head))
DyadSetToggleTemplate(UDyadSetToggle, StrictDyadSet, UTH(tail, head))

#endif // _ERGM_DYAD_HASHMAP_H_
