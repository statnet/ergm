/*  File inst/include/ergm_BDStratBlocks.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2022 Statnet Commons
 */
#ifndef _ERGM_BDSTRATBLOCKS_H_
#define _ERGM_BDSTRATBLOCKS_H_

#include "ergm_block.h"
#include "ergm_BDNodeLists.h"

typedef struct {
  Block ***blocks;
  int strat_nmixtypes;
  int *nblocks;

  BDNodeLists *lists;
} BDStratBlocks;

static inline BDStratBlocks *BDStratBlocksInitialize(BDNodeLists *lists,
                                                     int strat_nlevels,
                                                     int strat_nmixtypes,
                                                     int *strat_tails,
                                                     int *strat_heads,
                                                     int blocks_nlevels,
                                                     int *blocks_mixtypes,
                                                     int *blocks_tails,
                                                     int *blocks_heads,
                                                     int bd_nlevels,
                                                     int *bd_mixtypes,
                                                     int *bd_tails,
                                                     int *bd_heads,
                                                     Network *nwp) {
  BDStratBlocks *blocks = Calloc(1, BDStratBlocks);

  // do some copying
  blocks->lists = lists;

  // set up blocks
  blocks->blocks = Calloc(strat_nmixtypes, Block **);
  blocks->nblocks = Calloc(strat_nmixtypes, int);
  blocks->strat_nmixtypes = strat_nmixtypes;
  for(int i = 0; i < strat_nmixtypes; i++) {
    int strat_diag = (strat_tails[i] == strat_heads[i]);

    int nblocksmixtypes = blocks_mixtypes[strat_diag];

    int nblocksoffdiag = blocks_mixtypes[0] - blocks_mixtypes[1];
    int nblocksdiag = blocks_mixtypes[1] - nblocksoffdiag;

    int base_nblocks = (1 + !strat_diag)*nblocksoffdiag*bd_mixtypes[0] + nblocksdiag*bd_mixtypes[strat_diag];
    blocks->nblocks[i] = DIRECTED ? 4*base_nblocks : base_nblocks;
    blocks->blocks[i] = Calloc(blocks->nblocks[i], Block *);
    int l = 0;
    for(int j = 0; j < nblocksmixtypes; j++) {
      int blocks_diag = (blocks_tails[j] == blocks_heads[j]);

      for(int k = 0; k < bd_mixtypes[strat_diag && blocks_diag]; k++) {
        int tail_attr = strat_tails[i]*blocks_nlevels*bd_nlevels + blocks_tails[j]*bd_nlevels + bd_tails[k];
        int head_attr = strat_heads[i]*blocks_nlevels*bd_nlevels + blocks_heads[j]*bd_nlevels + bd_heads[k];
        if(DIRECTED) {
          blocks->blocks[i][l + 0] = BlockInitialize(lists->boths_by_attr[tail_attr][bd_heads[k]], lists->boths_by_attr[head_attr][bd_tails[k]], tail_attr == head_attr, DIRECTED);
          blocks->blocks[i][l + 1] = BlockInitialize(lists->tails_by_attr[tail_attr][bd_heads[k]], lists->boths_by_attr[head_attr][bd_tails[k]], FALSE, DIRECTED);
          blocks->blocks[i][l + 2] = BlockInitialize(lists->boths_by_attr[tail_attr][bd_heads[k]], lists->heads_by_attr[head_attr][bd_tails[k]], FALSE, DIRECTED);
          blocks->blocks[i][l + 3] = BlockInitialize(lists->tails_by_attr[tail_attr][bd_heads[k]], lists->heads_by_attr[head_attr][bd_tails[k]], FALSE, DIRECTED);
          l += 4;
        } else {
          blocks->blocks[i][l + 0] = BlockInitialize(lists->tails_by_attr[tail_attr][bd_heads[k]], lists->heads_by_attr[head_attr][bd_tails[k]], tail_attr == head_attr, DIRECTED);
          l += 1;
        }
      }
    }
  }

  return blocks;
}

static inline void BDStratBlocksDestroy(BDStratBlocks *blocks) {
  for(int i = 0; i < blocks->strat_nmixtypes; i++) {
    for(int j = 0; j < blocks->nblocks[i]; j++) {
      BlockDestroy(blocks->blocks[i][j]);
    }
    Free(blocks->blocks[i]);
  }
  Free(blocks->blocks);
  Free(blocks->nblocks);
  Free(blocks);
}

static inline Dyad BDStratBlocksDyadCount(BDStratBlocks *blocks, int stratmixingtype) {
  Dyad dyadcount = 0;
  for(int i = 0; i < blocks->nblocks[stratmixingtype]; i++) {
    dyadcount += BlockDyadCount(blocks->blocks[stratmixingtype][i]);
  }
  return dyadcount;
}

static inline int BDStratBlocksDyadCountPositive(BDStratBlocks *blocks, int stratmixingtype) {
  for(int i = 0; i < blocks->nblocks[stratmixingtype]; i++) {
    if(BlockDyadCountPositive(blocks->blocks[stratmixingtype][i])) {
      return TRUE;
    }
  }
  return FALSE;
}

static inline void BDStratBlocksGetRandWithCount(Vertex *tail, Vertex *head, BDStratBlocks *blocks, int stratmixingtype, Dyad dyadcount) {
  Dyad dyadindex = 2*dyadcount*unif_rand();

  for(int i = 0; i < blocks->nblocks[stratmixingtype]; i++) {
    Dyad dyadsthistype = 2*BlockDyadCount(blocks->blocks[stratmixingtype][i]);

    if(dyadindex < dyadsthistype) {
      BlockPut2Dyad(tail, head, dyadindex, blocks->blocks[stratmixingtype][i]);
      break;
    } else {
      dyadindex -= dyadsthistype;
    }
  }
}

static inline void BDStratBlocksGetRand(Vertex *tail, Vertex *head, BDStratBlocks *blocks, int stratmixingtype) {
  BDStratBlocksGetRandWithCount(tail, head, blocks, stratmixingtype, BDStratBlocksDyadCount(blocks, stratmixingtype));
}

static inline Dyad BDStratBlocksDyadCountOnToggle(Vertex tail, Vertex head, BDStratBlocks *blocks, int stratmixingtype, int change, int tailcondition, int headcondition) {
  BDNodeListsToggleIf(tail, head, blocks->lists, tailcondition, headcondition);
  Dyad dyadcount = BDStratBlocksDyadCount(blocks, stratmixingtype);
  BDNodeListsToggleIf(tail, head, blocks->lists, tailcondition, headcondition);
  return dyadcount;
}

#endif
