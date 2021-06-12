/*  File inst/include/ergm_BDStratBlocks.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2021 Statnet Commons
 */
#ifndef _ERGM_BDSTRATBLOCKS_H_
#define _ERGM_BDSTRATBLOCKS_H_

#include "ergm_block.h"

typedef struct {
  Block ***blocks;
  int strat_nmixtypes;
  int *nblocks;
  
  NodeList *****tails;
  NodeList *****heads;
  NodeList *****boths;
  
  int **tailpos;
  int **headpos;
  int **bothpos;
  
  int *strat_vattr;
  int *blocks_vattr;
  int *bd_vattr;
    
  int directed;

  int strat_nlevels;
  int blocks_nlevels;
  int bd_nlevels;
  
  NodeList *last_tails;
  NodeList *last_heads;
} BDStratBlocks;

static inline BDStratBlocks *BDStratBlocksInitialize(int **maxout, 
                                                     int **maxin, 
                                                     int *strat_vattr, 
                                                     int strat_nlevels, 
                                                     int strat_nmixtypes, 
                                                     int *strat_tails, 
                                                     int *strat_heads,
                                                     int *blocks_vattr, 
                                                     int blocks_nlevels, 
                                                     int *blocks_mixtypes, 
                                                     int *blocks_tails,
                                                     int *blocks_heads,
                                                     int *bd_vattr, 
                                                     int bd_nlevels, 
                                                     int *bd_mixtypes, 
                                                     int *bd_tails,
                                                     int *bd_heads,                                                     
                                                     int *jointattrcounts, 
                                                     int **indegree,
                                                     int **outdegree,
                                                     Network *nwp) {
  BDStratBlocks *blocks = Calloc(1, BDStratBlocks);

  // do some copying  
  blocks->directed = DIRECTED;

  blocks->strat_nlevels = strat_nlevels;
  blocks->strat_vattr = strat_vattr;  

  blocks->blocks_nlevels = blocks_nlevels;
  blocks->blocks_vattr = blocks_vattr;

  blocks->bd_nlevels = bd_nlevels;
  blocks->bd_vattr = bd_vattr;

  // decrement attribute pointers so node indices line up correctly
  blocks->strat_vattr--;
  blocks->blocks_vattr--;
  blocks->bd_vattr--;

  // set up node lists
  blocks->bothpos = Calloc(bd_nlevels, int *);
  blocks->tailpos = DIRECTED ? Calloc(bd_nlevels, int *) : blocks->bothpos;  
  blocks->headpos = DIRECTED ? Calloc(bd_nlevels, int *) : blocks->bothpos;  
  
  for(int i = 0; i < bd_nlevels; i++) {
    blocks->bothpos[i] = Calloc(N_NODES + 1, int);
    if(DIRECTED) {
      blocks->tailpos[i] = Calloc(N_NODES + 1, int);  
      blocks->headpos[i] = Calloc(N_NODES + 1, int); 
    }
  }
    
  blocks->boths = Calloc(strat_nlevels, NodeList ****);
  blocks->tails = DIRECTED ? Calloc(strat_nlevels, NodeList ****) : blocks->boths;
  blocks->heads = DIRECTED ? Calloc(strat_nlevels, NodeList ****) : blocks->boths;
  
  for(int i = 0; i < strat_nlevels; i++) {
    blocks->boths[i] = Calloc(blocks_nlevels, NodeList ***);
    if(DIRECTED) {
      blocks->tails[i] = Calloc(blocks_nlevels, NodeList ***);
      blocks->heads[i] = Calloc(blocks_nlevels, NodeList ***);
    }
    
    for(int j = 0; j < blocks_nlevels; j++) {
      blocks->boths[i][j] = Calloc(bd_nlevels, NodeList **);
      if(DIRECTED) {
        blocks->tails[i][j] = Calloc(bd_nlevels, NodeList **);
        blocks->heads[i][j] = Calloc(bd_nlevels, NodeList **);          
      }
      
      for(int k = 0; k < bd_nlevels; k++) {
        blocks->boths[i][j][k] = Calloc(bd_nlevels, NodeList *);
        if(DIRECTED) {
          blocks->tails[i][j][k] = Calloc(bd_nlevels, NodeList *);
          blocks->heads[i][j][k] = Calloc(bd_nlevels, NodeList *);          
        }
        
        for(int l = 0; l < bd_nlevels; l++) {
          blocks->boths[i][j][k][l] = NodeListInitialize(jointattrcounts[i*blocks_nlevels*bd_nlevels + j*bd_nlevels + k], blocks->bothpos[l]);
          if(DIRECTED) {
            blocks->tails[i][j][k][l] = NodeListInitialize(jointattrcounts[i*blocks_nlevels*bd_nlevels + j*bd_nlevels + k], blocks->tailpos[l]);
            blocks->heads[i][j][k][l] = NodeListInitialize(jointattrcounts[i*blocks_nlevels*bd_nlevels + j*bd_nlevels + k], blocks->headpos[l]);
          }
        }
      }
    }
  }

  for(Vertex vertex = 1; vertex <= N_NODES; vertex++) {
    int strat_val = blocks->strat_vattr[vertex];
    int blocks_val = blocks->blocks_vattr[vertex];
    int bd_val = blocks->bd_vattr[vertex];
    
    for(int i = 0; i < bd_nlevels; i++) {
      if(DIRECTED) {
        if(indegree[i][vertex] < maxin[i][vertex] && outdegree[i][vertex] < maxout[i][vertex]) {
          NodeListInsert(blocks->boths[strat_val][blocks_val][bd_val][i], vertex);        
        } else if(outdegree[i][vertex] < maxout[i][vertex]) {
          NodeListInsert(blocks->tails[strat_val][blocks_val][bd_val][i], vertex);
        } else if(indegree[i][vertex] < maxin[i][vertex]) {
          NodeListInsert(blocks->heads[strat_val][blocks_val][bd_val][i], vertex);
        }        
      } else if(indegree[i][vertex] + outdegree[i][vertex] < maxout[i][vertex]) {
        NodeListInsert(blocks->boths[strat_val][blocks_val][bd_val][i], vertex);        
      }
    }
  }

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
        if(DIRECTED) {
          blocks->blocks[i][l + 0] = BlockInitialize(blocks->boths[strat_tails[i]][blocks_tails[j]][bd_tails[k]][bd_heads[k]], blocks->boths[strat_heads[i]][blocks_heads[j]][bd_heads[k]][bd_tails[k]], bd_tails[k] == bd_heads[k] && blocks_diag && strat_diag, DIRECTED);
          blocks->blocks[i][l + 1] = BlockInitialize(blocks->tails[strat_tails[i]][blocks_tails[j]][bd_tails[k]][bd_heads[k]], blocks->boths[strat_heads[i]][blocks_heads[j]][bd_heads[k]][bd_tails[k]], FALSE, DIRECTED);
          blocks->blocks[i][l + 2] = BlockInitialize(blocks->boths[strat_tails[i]][blocks_tails[j]][bd_tails[k]][bd_heads[k]], blocks->heads[strat_heads[i]][blocks_heads[j]][bd_heads[k]][bd_tails[k]], FALSE, DIRECTED);
          blocks->blocks[i][l + 3] = BlockInitialize(blocks->tails[strat_tails[i]][blocks_tails[j]][bd_tails[k]][bd_heads[k]], blocks->heads[strat_heads[i]][blocks_heads[j]][bd_heads[k]][bd_tails[k]], FALSE, DIRECTED);          
          l += 4;
        } else {
          blocks->blocks[i][l + 0] = BlockInitialize(blocks->tails[strat_tails[i]][blocks_tails[j]][bd_tails[k]][bd_heads[k]], blocks->heads[strat_heads[i]][blocks_heads[j]][bd_heads[k]][bd_tails[k]], bd_tails[k] == bd_heads[k] && blocks_diag && strat_diag, DIRECTED);
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
  
  for(int i = 0; i < blocks->strat_nlevels; i++) {
    for(int j = 0; j < blocks->blocks_nlevels; j++) {
      for(int k = 0; k < blocks->bd_nlevels; k++) {
        for(int l = 0; l < blocks->bd_nlevels; l++) {
          NodeListDestroy(blocks->boths[i][j][k][l]);
          if(blocks->directed) {
            NodeListDestroy(blocks->tails[i][j][k][l]);
            NodeListDestroy(blocks->heads[i][j][k][l]);
          }            
        }
        Free(blocks->boths[i][j][k]);
        if(blocks->directed) {
          Free(blocks->tails[i][j][k]);
          Free(blocks->heads[i][j][k]);
        }                    
      }
      Free(blocks->boths[i][j]);
      if(blocks->directed) {
        Free(blocks->tails[i][j]);
        Free(blocks->heads[i][j]);
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
  
  for(int i = 0; i < blocks->bd_nlevels; i++) {
    Free(blocks->bothpos[i]);
    if(blocks->directed) {
      Free(blocks->tailpos[i]);
      Free(blocks->headpos[i]);
    }
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

static inline void BDStratBlocksSetLast(Vertex tail, Vertex head, int edgestate, BDStratBlocks *blocks) {
  if(edgestate) {
    if(blocks->directed && blocks->headpos[blocks->bd_vattr[head]][tail]) {
      blocks->last_tails = blocks->boths[blocks->strat_vattr[tail]][blocks->blocks_vattr[tail]][blocks->bd_vattr[tail]][blocks->bd_vattr[head]];
    } else {
      blocks->last_tails = blocks->tails[blocks->strat_vattr[tail]][blocks->blocks_vattr[tail]][blocks->bd_vattr[tail]][blocks->bd_vattr[head]];      
    }
    
    if(blocks->directed && blocks->tailpos[blocks->bd_vattr[tail]][head]) {
      blocks->last_heads = blocks->boths[blocks->strat_vattr[head]][blocks->blocks_vattr[head]][blocks->bd_vattr[head]][blocks->bd_vattr[tail]];
    } else {
      blocks->last_heads = blocks->heads[blocks->strat_vattr[head]][blocks->blocks_vattr[head]][blocks->bd_vattr[head]][blocks->bd_vattr[tail]];      
    }
  } else {
    if(blocks->directed && blocks->bothpos[blocks->bd_vattr[head]][tail]) {
      blocks->last_tails = blocks->boths[blocks->strat_vattr[tail]][blocks->blocks_vattr[tail]][blocks->bd_vattr[tail]][blocks->bd_vattr[head]];
    } else {
      blocks->last_tails = blocks->tails[blocks->strat_vattr[tail]][blocks->blocks_vattr[tail]][blocks->bd_vattr[tail]][blocks->bd_vattr[head]];      
    }
    
    if(blocks->directed && blocks->bothpos[blocks->bd_vattr[tail]][head]) {
      blocks->last_heads = blocks->boths[blocks->strat_vattr[head]][blocks->blocks_vattr[head]][blocks->bd_vattr[head]][blocks->bd_vattr[tail]];
    } else {
      blocks->last_heads = blocks->heads[blocks->strat_vattr[head]][blocks->blocks_vattr[head]][blocks->bd_vattr[head]][blocks->bd_vattr[tail]];      
    }      
  }
}

static inline void BDStratBlocksToggleIf(Vertex tail, Vertex head, BDStratBlocks *blocks, int tailcondition, int headcondition) {
  if(tailcondition) {
    if(blocks->directed && (blocks->bothpos[blocks->bd_vattr[head]][tail] || blocks->headpos[blocks->bd_vattr[head]][tail])) {
      NodeListToggle(blocks->heads[blocks->strat_vattr[tail]][blocks->blocks_vattr[tail]][blocks->bd_vattr[tail]][blocks->bd_vattr[head]], tail);
    }
    NodeListToggle(blocks->last_tails, tail);
  }

  if(headcondition) {
    if(blocks->directed && (blocks->bothpos[blocks->bd_vattr[tail]][head] || blocks->tailpos[blocks->bd_vattr[tail]][head])) {
      NodeListToggle(blocks->tails[blocks->strat_vattr[head]][blocks->blocks_vattr[head]][blocks->bd_vattr[head]][blocks->bd_vattr[tail]], head);
    }
    NodeListToggle(blocks->last_heads, head);
  }
}

static inline Dyad BDStratBlocksDyadCountOnToggle(Vertex tail, Vertex head, BDStratBlocks *blocks, int stratmixingtype, int change, int tailcondition, int headcondition) {
  if(tailcondition) {
    if(blocks->directed && (blocks->bothpos[blocks->bd_vattr[head]][tail] || blocks->headpos[blocks->bd_vattr[head]][tail])) {
      blocks->heads[blocks->strat_vattr[tail]][blocks->blocks_vattr[tail]][blocks->bd_vattr[tail]][blocks->bd_vattr[head]]->length -= change;
    }
    blocks->last_tails->length += change;
  }

  if(headcondition) {
    if(blocks->directed && (blocks->bothpos[blocks->bd_vattr[tail]][head] || blocks->tailpos[blocks->bd_vattr[tail]][head])) {
      blocks->tails[blocks->strat_vattr[head]][blocks->blocks_vattr[head]][blocks->bd_vattr[head]][blocks->bd_vattr[tail]]->length -= change;
    }
    blocks->last_heads->length += change;
  }

  Dyad dyadcount = BDStratBlocksDyadCount(blocks, stratmixingtype);

  if(tailcondition) {
    if(blocks->directed && (blocks->bothpos[blocks->bd_vattr[head]][tail] || blocks->headpos[blocks->bd_vattr[head]][tail])) {
      blocks->heads[blocks->strat_vattr[tail]][blocks->blocks_vattr[tail]][blocks->bd_vattr[tail]][blocks->bd_vattr[head]]->length += change;
    }
    blocks->last_tails->length -= change;
  }

  if(headcondition) {
    if(blocks->directed && (blocks->bothpos[blocks->bd_vattr[tail]][head] || blocks->tailpos[blocks->bd_vattr[tail]][head])) {
      blocks->tails[blocks->strat_vattr[head]][blocks->blocks_vattr[head]][blocks->bd_vattr[head]][blocks->bd_vattr[tail]]->length += change;
    }
    blocks->last_heads->length -= change;
  }

  return dyadcount;
}

#endif
