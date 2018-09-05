#ifndef _ERGM_DYAD_HASHMAP_H_
#define _ERGM_DYAD_HASHMAP_H_

#include "ergm_edgetree.h"
#include "ergm_khash.h"

/* Data structure to represent a dyad that can serve as a key to the hash. */
struct TailHead{
  Vertex tail, head;
};

/* Helper macros to construct TailHeads; the undirected variant ensures that t < h. */
#define TH(t,h,d) (d?THD(t,h):THU(t,h))
#define THD(t,h) ((struct TailHead){.tail=(t),.head=(h)})
#define THU(t,h) ((struct TailHead){.tail=MIN((t),(h)),.head=MAX((t),(h))})

/* Hash and comparison functions designed for tail-head pairs. */
// The of macro-ing this due to Bob Jenkins.
#define ROT_INT(x,k) (((x)<<(k)) | ((x)>>(32-(k))))
static inline unsigned int kh_scramble_int(unsigned int a){
  a ^= 0xc761c23c;
  a += ROT_INT(a,30);
  a ^= ROT_INT(a,11);
  a += ROT_INT(a,20);
  return a;
}
#define kh_vertexvertex_hash_func(key) (khint32_t)(kh_scramble_int(ROT_INT((key).tail,16) ^ (key).head))
#define kh_vertexvertex_hash_equal(a,b) (a.tail==b.tail && a.head==b.head)

/* Predefined khash type for mapping dyads onto unsigned ints. */
KHASH_INIT(DyadMapUInt, struct TailHead, unsigned int, TRUE, kh_vertexvertex_hash_func, kh_vertexvertex_hash_equal)
typedef khash_t(DyadMapUInt) StoreDyadMapUInt;

/* Utility function declarations. */
void PrintDyadMapUInt(StoreDyadMapUInt *h);

#endif // _ERGM_DYAD_HASHMAP_H_
