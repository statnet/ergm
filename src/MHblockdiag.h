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
} MH_BlockDiagSampInfo;

static inline MH_BlockDiagSampInfo unpack_BlockDiagSampInfo(double **inputs, Vertex bipartite, unsigned int directed_flag){
  double *x = *inputs;
  Vertex n = (x++)[0]; // number of blocks

  MH_BlockDiagSampInfo out = {
    .directed_flag=directed_flag,
    .n=n};

  out.epos=x; x+=n+1;

  if(bipartite){
    out.apos=x; x+=n+1;
  }else out.apos=out.epos;
  
  out.cwt=x; x+=n;

  *inputs=x;
  return out;
}

static inline void GetRandDyadBlockDiag(Vertex *tail, Vertex *head, const MH_BlockDiagSampInfo *b){
  Vertex blk = 1, t, h;
  double r = unif_rand();
  // FIXME: Change to a binary search to change from O(b->n) to O(log(b->n)).
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

static inline Vertex GetBlockID(Vertex tail, Vertex head, const MH_BlockDiagSampInfo *b){
  Vertex tblk = 1;
  while(tail>b->epos[tblk]) tblk++; 
  Vertex hblk = 1;
  while(head>b->epos[hblk]) hblk++;
  if(tblk==hblk) return tblk;
  else return 0;
}

#endif
