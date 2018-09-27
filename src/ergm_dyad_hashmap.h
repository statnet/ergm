#ifndef _ERGM_DYAD_HASHMAP_H_
#define _ERGM_DYAD_HASHMAP_H_

#include "stdbool.h"
#include "ergm_edgetree.h"

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

/* Helper macros to construct TailHeads; the undirected variant ensures that t < h. */
#define TH(t,h,d) (d?THD(t,h):THU(t,h))
#define THD(t,h) ((TailHead){.tail=(t),.head=(h)})
#define THU(t,h) ((TailHead){.tail=MIN((t),(h)),.head=MAX((t),(h))})

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
KHASH_INIT(DyadMapUInt, TailHead, unsigned int, true, kh_vertexvertex_hash_func, kh_vertexvertex_hash_equal)
typedef khash_t(DyadMapUInt) StoreDyadMapUInt;

/* Accessors, modifiers, and incrementors. */
#define _GETDMUI3(tail, head, hashmap) (kh_getval(DyadMapUInt, hashmap, TH(tail,head, DIRECTED), 0))
#define _GETDMUI4(tail, head, directed, hashmap) (kh_getval(DyadMapUInt, hashmap, TH(tail,head, directed), 0))
#define GETDMUI(...) _GET_OVERRIDE4(__VA_ARGS__, _GETDMUI4, _GETDMUI3,)(__VA_ARGS__)
#define SETDMUI4(tail, head, v, hashmap) {if(v==0) kh_unset(DyadMapUInt, hashmap, TH(tail,head, DIRECTED)); else kh_set(DyadMapUInt, hashmap, TH(tail,head, DIRECTED), v)}
#define SETDMUI5(tail, head, directed, v, hashmap) {if(v==0) kh_unset(DyadMapUInt, hashmap, TH(tail,head, directed)); else kh_set(DyadMapUInt, hashmap, TH(tail,head, directed), v)}
#define SETDMUI(...) _GET_OVERRIDE5(__VA_ARGS__, _SETDMUI5, _SETDMUI4,)(__VA_ARGS__)


static inline void IncDyadMapUInt(TailHead th, int inc, StoreDyadMapUInt *spcache){
  if(inc!=0){
    khiter_t pos = kh_get(DyadMapUInt, spcache, th);
    unsigned int val = pos==kh_none ? 0 : kh_value(spcache, pos);
    val += inc;
    if(val==0){
      kh_del(DyadMapUInt, spcache, pos);
    }else{
      if(pos==kh_none){
	int ret;
	pos = kh_put(DyadMapUInt, spcache, th, &ret);
      }
      kh_val(spcache, pos) = val;
    }
  }
}


/* Predefined khash type for dyad sets. This may or may not be faster than edgetree. */
KHASH_INIT(DyadSet, TailHead, char, false, kh_vertexvertex_hash_func, kh_vertexvertex_hash_equal)
typedef khash_t(DyadSet) StoreDyadSet;

// Toggle an element of a DyadSet.
static inline bool DyadSetToggle(TailHead th, StoreDyadSet *h){
  int ret;
  // Attempt insertion
  khiter_t i = kh_put(DyadSet, h, th, &ret);
  if(ret==0){
    // Already present: delete
    kh_del(DyadSet, h, i);
    return false;
  }else{
    // Inserted by kh_put above
    return true;
  }
}

/* Utility function declarations. */
void PrintDyadMapUInt(StoreDyadMapUInt *h);
void PrintDyadSet(StoreDyadSet *h);
StoreDyadSet *NetworkToDyadSet(Network *nwp);

#endif // _ERGM_DYAD_HASHMAP_H_
