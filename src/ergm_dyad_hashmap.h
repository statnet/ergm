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
#define kh_vertexvertex_hash_func(key) (khint32_t)(kh_scramble_int(ROT_INT((key).tail,16) ^ (key).head))
#define kh_vertexvertex_hash_equal(a,b) (a.tail==b.tail && a.head==b.head)

/* Predefined khash type for mapping dyads onto unsigned ints. */
KHASH_INIT(DyadMapUInt, struct TailHead, unsigned int, TRUE, kh_vertexvertex_hash_func, kh_vertexvertex_hash_equal)
typedef khash_t(DyadMapUInt) StoreDyadMapUInt;

/* Predefined khash type for dyad sets. This may or may not be faster than edgetree. */
KHASH_INIT(DyadSet, struct TailHead, char, FALSE, kh_vertexvertex_hash_func, kh_vertexvertex_hash_equal)
typedef khash_t(DyadSet) StoreDyadSet;

/* Utility function declarations. */
void PrintDyadMapUInt(StoreDyadMapUInt *h);

/* Utility function declarations. */
void PrintDyadSet(StoreDyadSet *h);

#endif // _ERGM_DYAD_HASHMAP_H_
