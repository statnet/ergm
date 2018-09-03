#ifndef _CHANGESTATS_SPCACHE_H_
#define _CHANGESTATS_SPCACHE_H_

#include "ergm_edgetree.h"
#include "khash_helpers.h"

struct TailHead{
  Vertex tail, head;
};

#define TH(t,h,d) (d?THD(t,h):THU(t,h))
#define THD(t,h) ((struct TailHead){.tail=(t),.head=(h)})
#define THU(t,h) ((struct TailHead){.tail=MIN((t),(h)),.head=MAX((t),(h))})

static inline unsigned int kh_scramble_int(unsigned int a){
  a ^= 0xc761c23c;
  a += a>>2 ^ a<<30;
  a ^= a>>21 ^ a<<11;
  a += a>>12 ^ a<<20;
  return a;
}

/* #define kh_vertexvertex_hash_func(key) kh_int64_hash_func((Dyad)key.tail << 32 | key.head) */
#define kh_vertexvertex_hash_func(key) (khint32_t)(kh_scramble_int((key).tail<<16 ^ (key).tail>>16 ^ (key).head))
#define kh_vertexvertex_hash_equal(a,b) (a.tail==b.tail && a.head==b.head)

KHASH_INIT(EdgeMapUInt, struct TailHead, unsigned int, TRUE, kh_vertexvertex_hash_func, kh_vertexvertex_hash_equal)
typedef khash_t(EdgeMapUInt) StoreEdgeMapUInt;

#endif
