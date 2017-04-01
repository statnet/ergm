/* Utilities to facilitate block-diagonal proposals. */

/* # Standard format for block-diagonal information for MHproposals #

   Let double *x point to the start of the block information segment
   and let b be the number of blocks. Here, the slice notation x[a :
   b] denotes elements of the array x[a] through x[b] inclusive.

   ## Unipartite networks ##

   x[0] = number of blocks (b)

   x[1 : 1+b] = a (b+1)-vector whose 0th element is 0 and whose i'th
     element is the index of the last vertex from block i (indexed
     from 1). That is, vertices x[1+i-1]+1 through x[1+i-1+1] belong
     to block i.

   x[2+b : b+2+b-1] = a b-vector of cumulative probabilities for
     selecting each block set in the R code. (That is, the 0th element
     is the probability of the first block, the last element is always
     1, etc.)

   ## Bipartite networks ##   

   x[0] = number of blocks (b)

   x[1 : 1+b] = a (b+1)-vector whose 0th element is 0 and whose i'th
     element is the index of the last b1 (ego/actor) vertex from block
     i (indexed from 1). That is, vertices x[1+i-1]+1 through
     x[1+i-1+1] belong to block i.

   x[1+b+1 : 1+2*b+1] = a (b+1)-vector whose 0th element is 0 and
     whose i'th element is the index of the last b2 (alter/event)
     vertex from block i (indexed from 1). That is, vertices
     x[1+b+1+i-1]+1 through x[1+b+1+i-1+1] belong to block i.

   x[3+2*b : 3+3*b-1] = a b-vector of cumulative probabilities for
     selecting each block set in the R code. (That is, the 0th element
     is the probability of the first block, the last element is always
     1, etc.)
*/

#ifndef MHBLOCKDIAG_H
#define MHBLOCKDIAG_H

#include "MHproposal.h"

typedef struct {
  double *epos; // starts and ends of blocks (b1)
  double *apos; // starts and ends of blocks (b2) == b1 if unipartite
  double *cwt; // cumulative selection probabilities
  Dyad ndyads; // number of dyads
  Vertex n; // number of blocks
  unsigned int directed_flag; // directed flag
} MH_BlockDiagInfo;

static inline MH_BlockDiagInfo unpack_BlockDiagInfo(double *inputs, Vertex bipartite, unsigned int directed_flag){
  MH_BlockDiagInfo out = {
    .directed_flag=directed_flag,
    .n=inputs[0],
    .epos=inputs+1,
    .apos=bipartite ? inputs+2+(Vertex)(inputs[1]) : inputs+1,
    .cwt=bipartite ? inputs+3+(Vertex)(inputs[1])+(Vertex)(inputs[1]) : inputs+2+(Vertex)(inputs[1])
  };
  return out;
}

static inline void GetRandDyadBlockDiag(Vertex *tail, Vertex *head, const MH_BlockDiagInfo *b){
  Vertex blk = 1, t, h;
  double r = unif_rand();
  while(r>b->cwt[blk-1]) blk++;
  t = b->epos[blk-1]+1 + unif_rand() * (b->epos[blk]-b->epos[blk-1]);
  while ((h = b->apos[blk-1]+1 + unif_rand() * (b->apos[blk]-b->apos[blk-1])) == t);
  if (!b->directed_flag && t > h) {
    *tail = h;
    *head = t;
  }else{
    *tail = t;
    *head = h;
  }
}

#endif
