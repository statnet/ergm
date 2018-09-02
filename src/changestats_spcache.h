#ifndef _CHANGESTATS_SPCACHE_H_
#define _CHANGESTATS_SPCACHE_H_

#include "ergm_edgetree.h"
#include "khash_helpers.h"

struct TailHead{
  Vertex tail, head;
};

#define TH(t,h,d) (d?(struct TailHead){.tail=(t),.head=(h)}:(struct TailHead){.tail=MIN((t),(h)),.head=MAX((t),(h))})

#define kh_vertexvertex_hash_func(key) kh_int64_hash_func((Dyad)key.tail << 32 | key.head)
#define kh_vertexvertex_hash_equal(a,b) (a.tail==b.tail && a.head==b.head)

KHASH_INIT(EdgeMapUInt, struct TailHead, unsigned int, TRUE, kh_vertexvertex_hash_func, kh_vertexvertex_hash_equal)
typedef khash_t(EdgeMapUInt) StoreEdgeMapUInt;

#endif
