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

typedef struct {
  Block ***blocks;
  int strat_nmixtypes;
  int *nblocks;
  
  NodeList *****tails_by_attr;
  NodeList *****heads_by_attr;
  NodeList *****boths_by_attr;
  
  NodeList ***tails;
  NodeList ***heads;
  NodeList ***boths;

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
  
  int **maxout;
  int **maxin;
  int **minout;
  int **minin;

  int **indegree;
  int **outdegree;
  
  Network *nwp;
} BDStratBlocks;

static inline BDStratBlocks *BDStratBlocksInitialize(int **maxout, 
                                                     int **maxin, 
                                                     int **minout, 
                                                     int **minin,
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
  blocks->nwp = nwp;
  
  blocks->maxout = maxout;
  blocks->maxin = maxin;
  blocks->minout = minout;
  blocks->minin = minin;

  blocks->indegree = indegree;
  blocks->outdegree = outdegree;
  
  blocks->directed = DIRECTED;
  
  blocks->strat_nlevels = strat_nlevels;
  blocks->strat_vattr = strat_vattr;  

  blocks->blocks_nlevels = blocks_nlevels;
  blocks->blocks_vattr = blocks_vattr;

  blocks->bd_nlevels = bd_nlevels;
  blocks->bd_vattr = bd_vattr;

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
    
  blocks->boths_by_attr = Calloc(strat_nlevels, NodeList ****);
  blocks->tails_by_attr = DIRECTED ? Calloc(strat_nlevels, NodeList ****) : blocks->boths_by_attr;
  blocks->heads_by_attr = DIRECTED ? Calloc(strat_nlevels, NodeList ****) : blocks->boths_by_attr;
  
  for(int i = 0; i < strat_nlevels; i++) {
    blocks->boths_by_attr[i] = Calloc(blocks_nlevels, NodeList ***);
    if(DIRECTED) {
      blocks->tails_by_attr[i] = Calloc(blocks_nlevels, NodeList ***);
      blocks->heads_by_attr[i] = Calloc(blocks_nlevels, NodeList ***);
    }
    
    for(int j = 0; j < blocks_nlevels; j++) {
      blocks->boths_by_attr[i][j] = Calloc(bd_nlevels, NodeList **);
      if(DIRECTED) {
        blocks->tails_by_attr[i][j] = Calloc(bd_nlevels, NodeList **);
        blocks->heads_by_attr[i][j] = Calloc(bd_nlevels, NodeList **);          
      }
      
      for(int k = 0; k < bd_nlevels; k++) {
        blocks->boths_by_attr[i][j][k] = Calloc(bd_nlevels, NodeList *);
        if(DIRECTED) {
          blocks->tails_by_attr[i][j][k] = Calloc(bd_nlevels, NodeList *);
          blocks->heads_by_attr[i][j][k] = Calloc(bd_nlevels, NodeList *);          
        }
        
        for(int l = 0; l < bd_nlevels; l++) {
          blocks->boths_by_attr[i][j][k][l] = NodeListInitialize(jointattrcounts[i*blocks_nlevels*bd_nlevels + j*bd_nlevels + k], blocks->bothpos[l]);
          if(DIRECTED) {
            blocks->tails_by_attr[i][j][k][l] = NodeListInitialize(jointattrcounts[i*blocks_nlevels*bd_nlevels + j*bd_nlevels + k], blocks->tailpos[l]);
            blocks->heads_by_attr[i][j][k][l] = NodeListInitialize(jointattrcounts[i*blocks_nlevels*bd_nlevels + j*bd_nlevels + k], blocks->headpos[l]);
          }
        }
      }
    }
  }

  blocks->boths = Calloc(N_NODES + 1, NodeList **);
  blocks->tails = DIRECTED ? Calloc(N_NODES + 1, NodeList **) : blocks->boths;
  blocks->heads = DIRECTED ? Calloc(N_NODES + 1, NodeList **) : blocks->boths;

  for(Vertex vertex = 1; vertex <= N_NODES; vertex++) {
    int strat_val = blocks->strat_vattr[vertex];
    int blocks_val = blocks->blocks_vattr[vertex];
    int bd_val = blocks->bd_vattr[vertex];
    
    for(int i = 0; i < bd_nlevels; i++) {
      if(DIRECTED) {
        if(indegree[i][vertex] < maxin[i][vertex] && outdegree[i][vertex] < maxout[i][vertex]) {
          NodeListInsert(blocks->boths_by_attr[strat_val][blocks_val][bd_val][i], vertex);        
        } else if(outdegree[i][vertex] < maxout[i][vertex]) {
          NodeListInsert(blocks->tails_by_attr[strat_val][blocks_val][bd_val][i], vertex);
        } else if(indegree[i][vertex] < maxin[i][vertex]) {
          NodeListInsert(blocks->heads_by_attr[strat_val][blocks_val][bd_val][i], vertex);
        }        
      } else if(indegree[i][vertex] + outdegree[i][vertex] < maxout[i][vertex]) {
        NodeListInsert(blocks->boths_by_attr[strat_val][blocks_val][bd_val][i], vertex);        
      }
    }
    
    blocks->boths[vertex] = blocks->boths_by_attr[strat_val][blocks_val][bd_val];
    if(DIRECTED) {
      blocks->tails[vertex] = blocks->tails_by_attr[strat_val][blocks_val][bd_val];
      blocks->heads[vertex] = blocks->heads_by_attr[strat_val][blocks_val][bd_val];                
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
          blocks->blocks[i][l + 0] = BlockInitialize(blocks->boths_by_attr[strat_tails[i]][blocks_tails[j]][bd_tails[k]][bd_heads[k]], blocks->boths_by_attr[strat_heads[i]][blocks_heads[j]][bd_heads[k]][bd_tails[k]], bd_tails[k] == bd_heads[k] && blocks_diag && strat_diag, DIRECTED);
          blocks->blocks[i][l + 1] = BlockInitialize(blocks->tails_by_attr[strat_tails[i]][blocks_tails[j]][bd_tails[k]][bd_heads[k]], blocks->boths_by_attr[strat_heads[i]][blocks_heads[j]][bd_heads[k]][bd_tails[k]], FALSE, DIRECTED);
          blocks->blocks[i][l + 2] = BlockInitialize(blocks->boths_by_attr[strat_tails[i]][blocks_tails[j]][bd_tails[k]][bd_heads[k]], blocks->heads_by_attr[strat_heads[i]][blocks_heads[j]][bd_heads[k]][bd_tails[k]], FALSE, DIRECTED);
          blocks->blocks[i][l + 3] = BlockInitialize(blocks->tails_by_attr[strat_tails[i]][blocks_tails[j]][bd_tails[k]][bd_heads[k]], blocks->heads_by_attr[strat_heads[i]][blocks_heads[j]][bd_heads[k]][bd_tails[k]], FALSE, DIRECTED);          
          l += 4;
        } else {
          blocks->blocks[i][l + 0] = BlockInitialize(blocks->tails_by_attr[strat_tails[i]][blocks_tails[j]][bd_tails[k]][bd_heads[k]], blocks->heads_by_attr[strat_heads[i]][blocks_heads[j]][bd_heads[k]][bd_tails[k]], bd_tails[k] == bd_heads[k] && blocks_diag && strat_diag, DIRECTED);
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
          NodeListDestroy(blocks->boths_by_attr[i][j][k][l]);
          if(blocks->directed) {
            NodeListDestroy(blocks->tails_by_attr[i][j][k][l]);
            NodeListDestroy(blocks->heads_by_attr[i][j][k][l]);
          }            
        }
        Free(blocks->boths_by_attr[i][j][k]);
        if(blocks->directed) {
          Free(blocks->tails_by_attr[i][j][k]);
          Free(blocks->heads_by_attr[i][j][k]);
        }                    
      }
      Free(blocks->boths_by_attr[i][j]);
      if(blocks->directed) {
        Free(blocks->tails_by_attr[i][j]);
        Free(blocks->heads_by_attr[i][j]);
      }                    
    }
    Free(blocks->boths_by_attr[i]);
    if(blocks->directed) {
      Free(blocks->tails_by_attr[i]);
      Free(blocks->heads_by_attr[i]);
    }                    
  }
  Free(blocks->boths_by_attr);
  if(blocks->directed) {
    Free(blocks->tails_by_attr);
    Free(blocks->heads_by_attr);
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
  do {
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
  } while ((blocks->minout || blocks->minin) && 
           EdgetreeSearch(*tail, *head, blocks->nwp->outedges) &&
           ((blocks->minout && (blocks->minout[blocks->bd_vattr[*head]][*tail] == (blocks->nwp->directed_flag ? blocks->outdegree[blocks->bd_vattr[*head]][*tail] : blocks->outdegree[blocks->bd_vattr[*head]][*tail] + blocks->indegree[blocks->bd_vattr[*head]][*tail]))) ||
            (blocks->minin && (blocks->minin[blocks->bd_vattr[*tail]][*head] == (blocks->nwp->directed_flag ? blocks->indegree[blocks->bd_vattr[*tail]][*head] : blocks->outdegree[blocks->bd_vattr[*tail]][*head] + blocks->indegree[blocks->bd_vattr[*tail]][*head])))));  
}

static inline void BDStratBlocksGetRand(Vertex *tail, Vertex *head, BDStratBlocks *blocks, int stratmixingtype) {
  BDStratBlocksGetRandWithCount(tail, head, blocks, stratmixingtype, BDStratBlocksDyadCount(blocks, stratmixingtype));
}

static inline void BDStratBlocksSetLast(Vertex tail, Vertex head, int edgestate, BDStratBlocks *blocks) {
  if(edgestate) {
    if(blocks->directed && blocks->headpos[blocks->bd_vattr[head]][tail]) {
      blocks->last_tails = blocks->boths[tail][blocks->bd_vattr[head]];
    } else {
      blocks->last_tails = blocks->tails[tail][blocks->bd_vattr[head]];      
    }
    
    if(blocks->directed && blocks->tailpos[blocks->bd_vattr[tail]][head]) {
      blocks->last_heads = blocks->boths[head][blocks->bd_vattr[tail]];
    } else {
      blocks->last_heads = blocks->heads[head][blocks->bd_vattr[tail]];      
    }
  } else {
    if(blocks->directed && blocks->bothpos[blocks->bd_vattr[head]][tail]) {
      blocks->last_tails = blocks->boths[tail][blocks->bd_vattr[head]];
    } else {
      blocks->last_tails = blocks->tails[tail][blocks->bd_vattr[head]];      
    }
    
    if(blocks->directed && blocks->bothpos[blocks->bd_vattr[tail]][head]) {
      blocks->last_heads = blocks->boths[head][blocks->bd_vattr[tail]];
    } else {
      blocks->last_heads = blocks->heads[head][blocks->bd_vattr[tail]];      
    }      
  }
}

static inline void BDStratBlocksToggleIf(Vertex tail, Vertex head, BDStratBlocks *blocks, int tailcondition, int headcondition) {
  if(tailcondition) {
    if(blocks->directed && (blocks->bothpos[blocks->bd_vattr[head]][tail] || blocks->headpos[blocks->bd_vattr[head]][tail])) {
      NodeListToggle(blocks->heads[tail][blocks->bd_vattr[head]], tail);
    }
    NodeListToggle(blocks->last_tails, tail);
  }

  if(headcondition) {
    if(blocks->directed && (blocks->bothpos[blocks->bd_vattr[tail]][head] || blocks->tailpos[blocks->bd_vattr[tail]][head])) {
      NodeListToggle(blocks->tails[head][blocks->bd_vattr[tail]], head);
    }
    NodeListToggle(blocks->last_heads, head);
  }
}

static inline Dyad BDStratBlocksDyadCountOnToggle(Vertex tail, Vertex head, BDStratBlocks *blocks, int stratmixingtype, int change, int tailcondition, int headcondition) {
  if(tailcondition) {
    if(blocks->directed && (blocks->bothpos[blocks->bd_vattr[head]][tail] || blocks->headpos[blocks->bd_vattr[head]][tail])) {
      blocks->heads[tail][blocks->bd_vattr[head]]->length -= change;
    }
    blocks->last_tails->length += change;
  }

  if(headcondition) {
    if(blocks->directed && (blocks->bothpos[blocks->bd_vattr[tail]][head] || blocks->tailpos[blocks->bd_vattr[tail]][head])) {
      blocks->tails[head][blocks->bd_vattr[tail]]->length -= change;
    }
    blocks->last_heads->length += change;
  }

  Dyad dyadcount = BDStratBlocksDyadCount(blocks, stratmixingtype);

  if(tailcondition) {
    if(blocks->directed && (blocks->bothpos[blocks->bd_vattr[head]][tail] || blocks->headpos[blocks->bd_vattr[head]][tail])) {
      blocks->heads[tail][blocks->bd_vattr[head]]->length += change;
    }
    blocks->last_tails->length -= change;
  }

  if(headcondition) {
    if(blocks->directed && (blocks->bothpos[blocks->bd_vattr[tail]][head] || blocks->tailpos[blocks->bd_vattr[tail]][head])) {
      blocks->tails[head][blocks->bd_vattr[tail]]->length += change;
    }
    blocks->last_heads->length -= change;
  }

  return dyadcount;
}

#endif
