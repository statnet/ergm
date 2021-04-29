#ifndef _ERGM_BDSTRATBLOCKS_H_
#define _ERGM_BDSTRATBLOCKS_H_

#include "ergm_block.h"

typedef struct {
  Block ***blocks;
  int strat_nmixtypes;
  int *nblocks;
  
  NodeList ***tails;
  NodeList ***heads;
  NodeList ***boths;
  
  int *tailpos;
  int *headpos;
  int *bothpos;
  
  int *strat_vattr;
  int *bd_vattr;
    
  int directed;

  int strat_nlevels;
  int bd_nlevels;
} BDStratBlocks;

static inline BDStratBlocks *BDStratBlocksInitialize(int maxout, 
                                                     int maxin, 
                                                     int *strat_vattr, 
                                                     int strat_nlevels, 
                                                     int strat_nmixtypes, 
                                                     int *strat_tails, 
                                                     int *strat_heads,
                                                     int *bd_vattr, 
                                                     int bd_nlevels, 
                                                     int *bd_mixtypes, 
                                                     int *bd_tails,
                                                     int *bd_heads,
                                                     int *jointattrcounts, 
                                                     Network *nwp) {
  BDStratBlocks *blocks = Calloc(1, BDStratBlocks);

  // do some copying  
  blocks->directed = DIRECTED;

  blocks->strat_nlevels = strat_nlevels;
  blocks->strat_vattr = strat_vattr;  

  blocks->bd_nlevels = bd_nlevels;
  blocks->bd_vattr = bd_vattr;

  // decrement attribute pointers so node indices line up correctly
  blocks->strat_vattr--;
  blocks->bd_vattr--;

  // set up node lists
  blocks->bothpos = Calloc(N_NODES + 1, int);
  blocks->tailpos = DIRECTED ? Calloc(N_NODES + 1, int) : blocks->bothpos;  
  blocks->headpos = DIRECTED ? Calloc(N_NODES + 1, int) : blocks->bothpos;  
  
  blocks->boths = Calloc(strat_nlevels, NodeList **);
  blocks->tails = DIRECTED ? Calloc(strat_nlevels, NodeList **) : blocks->boths;
  blocks->heads = DIRECTED ? Calloc(strat_nlevels, NodeList **) : blocks->boths;
  
  for(int i = 0; i < strat_nlevels; i++) {
    blocks->boths[i] = Calloc(bd_nlevels, NodeList *);
    if(DIRECTED) {
      blocks->tails[i] = Calloc(bd_nlevels, NodeList *);
      blocks->heads[i] = Calloc(bd_nlevels, NodeList *);
    }
    for(int j = 0; j < bd_nlevels; j++) {
      blocks->boths[i][j] = NodeListInitialize(jointattrcounts[i*bd_nlevels + j], blocks->bothpos);
      if(DIRECTED) {
        blocks->tails[i][j] = NodeListInitialize(jointattrcounts[i*bd_nlevels + j], blocks->tailpos);
        blocks->heads[i][j] = NodeListInitialize(jointattrcounts[i*bd_nlevels + j], blocks->headpos);
      }
    }
  }
  
  for(Vertex vertex = 1; vertex <= N_NODES; vertex++) {
    int strat_val = blocks->strat_vattr[vertex];
    int bd_val = blocks->bd_vattr[vertex];
    
    if(DIRECTED) {
      if(IN_DEG[vertex] < maxin && OUT_DEG[vertex] < maxout) {
        NodeListInsert(blocks->boths[strat_val][bd_val], vertex);        
      } else if(OUT_DEG[vertex] < maxout) {
        NodeListInsert(blocks->tails[strat_val][bd_val], vertex);
      } else if(IN_DEG[vertex] < maxin) {
        NodeListInsert(blocks->heads[strat_val][bd_val], vertex);
      }        
    } else if(IN_DEG[vertex] + OUT_DEG[vertex] < maxout) {
      NodeListInsert(blocks->boths[strat_val][bd_val], vertex);        
    }
  }
  
  // set up blocks
  blocks->blocks = Calloc(strat_nmixtypes, Block **);
  blocks->nblocks = Calloc(strat_nmixtypes, int);
  blocks->strat_nmixtypes = strat_nmixtypes;
  for(int i = 0; i < strat_nmixtypes; i++) {
    int strat_diag = (strat_tails[i] == strat_heads[i]);
    blocks->nblocks[i] = DIRECTED ? 4*bd_mixtypes[strat_diag] : bd_mixtypes[strat_diag];
    blocks->blocks[i] = Calloc(blocks->nblocks[i], Block *);
    for(int j = 0; j < bd_mixtypes[strat_diag]; j++) {
      if(DIRECTED) {
        blocks->blocks[i][4*j + 0] = BlockInitialize(blocks->boths[strat_tails[i]][bd_tails[j]], blocks->boths[strat_heads[i]][bd_heads[j]], (bd_tails[j] == bd_heads[j]) && strat_diag, DIRECTED);
        blocks->blocks[i][4*j + 1] = BlockInitialize(blocks->tails[strat_tails[i]][bd_tails[j]], blocks->boths[strat_heads[i]][bd_heads[j]], FALSE, DIRECTED);
        blocks->blocks[i][4*j + 2] = BlockInitialize(blocks->boths[strat_tails[i]][bd_tails[j]], blocks->heads[strat_heads[i]][bd_heads[j]], FALSE, DIRECTED);
        blocks->blocks[i][4*j + 3] = BlockInitialize(blocks->tails[strat_tails[i]][bd_tails[j]], blocks->heads[strat_heads[i]][bd_heads[j]], FALSE, DIRECTED);          
      } else {
        blocks->blocks[i][j] = BlockInitialize(blocks->tails[strat_tails[i]][bd_tails[j]], blocks->heads[strat_heads[i]][bd_heads[j]], (bd_tails[j] == bd_heads[j]) && strat_diag, DIRECTED);
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
  
  for(int i = 0; i < blocks->strat_nlevels; i++) {
    for(int j = 0; j < blocks->bd_nlevels; j++) {
      NodeListDestroy(blocks->boths[i][j]);
      if(blocks->directed) {
        NodeListDestroy(blocks->tails[i][j]);
        NodeListDestroy(blocks->heads[i][j]);
      }
    }
    Free(blocks->boths[i]);
    if(blocks->directed) {
      Free(blocks->tails[i]);
      Free(blocks->heads[i]);
    }
  }
  Free(blocks->boths);
  if(blocks->directed) {
    Free(blocks->tails);
    Free(blocks->heads);
  }
  
  Free(blocks->bothpos);
  if(blocks->directed) {
    Free(blocks->tailpos);
    Free(blocks->headpos);
  }
    
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
      return;
    } else {
      dyadindex -= dyadsthistype;  
    }
  }
}

static inline void BDStratBlocksGetRand(Vertex *tail, Vertex *head, BDStratBlocks *blocks, int stratmixingtype) {
  BDStratBlocksGetRandWithCount(tail, head, blocks, stratmixingtype, BDStratBlocksDyadCount(blocks, stratmixingtype));
}

static inline void BDStratBlocksToggleIf(Vertex tail, Vertex head, BDStratBlocks *blocks, int tailcondition, int headcondition) {
  if(tailcondition) {
    if(blocks->directed && (blocks->bothpos[tail] || blocks->headpos[tail])) {
      NodeListToggle(blocks->boths[blocks->strat_vattr[tail]][blocks->bd_vattr[tail]], tail);
      NodeListToggle(blocks->heads[blocks->strat_vattr[tail]][blocks->bd_vattr[tail]], tail);
    } else {
      NodeListToggle(blocks->tails[blocks->strat_vattr[tail]][blocks->bd_vattr[tail]], tail);
    }
  }

  if(headcondition) {
    if(blocks->directed && (blocks->bothpos[head] || blocks->tailpos[head])) {
      NodeListToggle(blocks->boths[blocks->strat_vattr[head]][blocks->bd_vattr[head]], head);
      NodeListToggle(blocks->tails[blocks->strat_vattr[head]][blocks->bd_vattr[head]], head);
    } else {
      NodeListToggle(blocks->heads[blocks->strat_vattr[head]][blocks->bd_vattr[head]], head);
    }
  }
}

static inline Dyad BDStratBlocksDyadCountOnToggle(Vertex tail, Vertex head, BDStratBlocks *blocks, int stratmixingtype, int tailcondition, int headcondition) {
  BDStratBlocksToggleIf(tail, head, blocks, tailcondition, headcondition);
  Dyad dyadcount = BDStratBlocksDyadCount(blocks, stratmixingtype);
  BDStratBlocksToggleIf(tail, head, blocks, tailcondition, headcondition);
  return dyadcount;
}

#endif
