/*  File inst/include/ergm_block.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2022 Statnet Commons
 */
#ifndef _ERGM_BLOCK_H_
#define _ERGM_BLOCK_H_

#include "ergm_nodelist.h"

/*
   This data structure represents a block whose margins are NodeLists, which may
   be dynamic. It supports calculation of dyad count (including directedness and
   diagonal information) and extraction of a dyad at a specified doubled index.
   (Indices are doubled to keep things rectangular for diagonal, undirected blocks,
   whose unique dyads form a triangle.)
*/


typedef struct {
  NodeList *tails;
  NodeList *heads;
  int diagonal;
  int directed;
} Block;

static inline Block *BlockInitialize(NodeList *tails, NodeList *heads, int diagonal, int directed) {
  Block *block = Calloc(1, Block);
  block->tails = tails;
  block->heads = heads;
  block->diagonal = diagonal;
  block->directed = directed;
  return block;
}

static inline void BlockDestroy(Block *block) {
  Free(block);
}

static inline void BlockPut2Dyad(Vertex *tail, Vertex *head, Dyad dyadindex, Block *block) {
  int tailindex;
  int headindex;

  if(block->diagonal) {
    if(block->directed) {
      dyadindex /= 2;
    }
    tailindex = dyadindex / block->tails->length;
    headindex = dyadindex % (block->heads->length - 1);
    if(tailindex == headindex) {
      headindex = block->heads->length - 1;
    }
  } else {
    dyadindex /= 2;
    tailindex = dyadindex / block->heads->length;
    headindex = dyadindex % block->heads->length;
  }

  // 1-based indexing in NodeLists
  tailindex++;
  headindex++;

  if(block->tails->nodes[tailindex] < block->heads->nodes[headindex] || block->directed) {
    *tail = block->tails->nodes[tailindex];
    *head = block->heads->nodes[headindex];
  } else {
    *tail = block->heads->nodes[headindex];
    *head = block->tails->nodes[tailindex];
  }
}

static inline Dyad BlockDyadCount(Block *block) {
  if(block->diagonal) {
    if(block->directed) {
      return (Dyad)block->tails->length*(block->heads->length - 1);
    } else {
      return (Dyad)block->tails->length*(block->heads->length - 1)/2;
    }
  } else {
    return (Dyad)block->tails->length*block->heads->length;
  }
}

static inline int BlockDyadCountPositive(Block *block) {
  if(block->diagonal) {
    return block->tails->length > 1;
  } else {
    return block->tails->length > 0 && block->heads->length > 0;
  }
}

#endif
