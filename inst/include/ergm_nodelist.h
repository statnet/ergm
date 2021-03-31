#ifndef _ERGM_NODELIST_H_
#define _ERGM_NODELIST_H_

#include "ergm_MHproposal.h"
#include "ergm_edgetree_types.h"
#include "ergm_changestat.h"

typedef struct {
  Vertex *nodes;
  int length;
  int *nodepos;
} NL;

static inline NL *NLInitialize(int nnodes, int *nodepos) {
  NL *nodelist = Calloc(1, NL);
  nodelist->nodes = Calloc(nnodes + 1, Vertex);
  nodelist->nodepos = nodepos;
  return nodelist;
}

static inline void NLDestroy(NL *nodelist) {
  Free(nodelist->nodes);
  Free(nodelist);
}

static inline void NLInsert(NL *nodelist, Vertex node) {
  nodelist->length++;
  nodelist->nodes[nodelist->length] = node;
  nodelist->nodepos[node] = nodelist->length;  
}

static inline void NLDelete(NL *nodelist, Vertex node) {
  nodelist->nodes[nodelist->nodepos[node]] = nodelist->nodes[nodelist->length];
  nodelist->nodepos[nodelist->nodes[nodelist->length]] = nodelist->nodepos[node];
  nodelist->nodepos[node] = 0;
  nodelist->length--;
}

static inline void NLToggleKnown(NL *nodelist, Vertex node, int nodeflag) {
  if(nodeflag) {
    NLDelete(nodelist, node);
  } else {
    NLInsert(nodelist, node);
  }
}

static inline void NLPut2Dyad(Vertex *tail, Vertex *head, NL *tails, NL *heads, Dyad dyadindex, int diagonal, int directed) {
  int tailcount = tails->length;
  int headcount = heads->length;

  int tailindex;
  int headindex;
  
  if(diagonal) {
    if(directed) {
      dyadindex /= 2;
    }
    tailindex = dyadindex / tailcount;
    headindex = dyadindex % (headcount - 1);
    if(tailindex == headindex) {
      headindex = headcount - 1;
    }                  
  } else {
    dyadindex /= 2;
    tailindex = dyadindex / headcount;
    headindex = dyadindex % headcount;        
  }
  
  // 1-based indexing in NLs
  tailindex++;
  headindex++;
  
  if(tails->nodes[tailindex] < heads->nodes[headindex] || directed) {
    *tail = tails->nodes[tailindex];
    *head = heads->nodes[headindex];
  } else {
    *tail = heads->nodes[headindex];
    *head = tails->nodes[tailindex];      
  }
}

static inline Dyad NL2DyadCount(NL *one, NL *two, int diagonal, int directed) {
  if(diagonal) {
    if(directed) {
      return (Dyad)2*one->length*(two->length - 1);
    } else {
      return (Dyad)one->length*(two->length - 1);      
    }
  } else {
    return (Dyad)2*one->length*two->length;
  }
}

static inline Dyad NLDyadCount(NL *one, NL *two, int diagonal, int directed) {
  if(diagonal) {
    if(directed) {
      return (Dyad)one->length*(two->length - 1);
    } else {
      return (Dyad)one->length*(two->length - 1)/2;      
    }
  } else {
    return (Dyad)one->length*two->length;
  }
}

static inline int NLDyadCountPositive(NL *one, NL *two, int diagonal) {
  if(diagonal) {
    return one->length > 1;
  } else {
    return one->length > 0 && two->length > 0;
  }
}

typedef struct {
  NL ***tails;
  NL ***heads;
  NL ***boths;
  
  int *tailpos;
  int *headpos;
  int *bothpos;
  
  int *strat_tails;
  int *strat_heads;
  
  int *bd_tails;
  int *bd_heads;

  int *strat_vattr;
  int *bd_vattr;
  
  int *lengths;
  
  int *strat_diagonal;
  int *bd_diagonal;
  
  int directed;

  int strat_nlevels;
  int bd_nlevels;
} NodeList;

static inline NodeList *NodeListInitialize(int maxout, int maxin, int *strat_vattr, int strat_nlevels, int strat_nmixtypes, int *strat_tails, int *strat_heads,
                                           int *bd_vattr, int bd_nlevels, int *bd_mixtypes, int *bd_tails, int *bd_heads,
                                           int *jointattrcounts, Network *nwp) {
  NodeList *nodelist = Calloc(1, NodeList);
  nodelist->directed = DIRECTED;

  nodelist->strat_diagonal = Calloc(strat_nmixtypes, int);
  for(int i = 0; i < strat_nmixtypes; i++) {
    nodelist->strat_diagonal[i] = (strat_tails[i] == strat_heads[i]);
  }

  nodelist->bd_diagonal = Calloc(bd_mixtypes[0], int);
  for(int i = 0; i < bd_mixtypes[0]; i++) {
    nodelist->bd_diagonal[i] = (bd_tails[i] == bd_heads[i]);
  }

  nodelist->strat_tails = strat_tails;
  nodelist->strat_heads = strat_heads;

  nodelist->bd_tails = bd_tails;
  nodelist->bd_heads = bd_heads;

  nodelist->lengths = bd_mixtypes;

  nodelist->strat_vattr = strat_vattr - 1; // decrement so node indices line up correctly  
  nodelist->bd_vattr = bd_vattr - 1; // decrement so node indices line up correctly

  nodelist->strat_nlevels = strat_nlevels;
  nodelist->bd_nlevels = bd_nlevels;

  nodelist->bothpos = Calloc(N_NODES + 1, int);
  nodelist->tailpos = DIRECTED ? Calloc(N_NODES + 1, int) : nodelist->bothpos;  
  nodelist->headpos = DIRECTED ? Calloc(N_NODES + 1, int) : nodelist->bothpos;  
  
  nodelist->boths = Calloc(strat_nlevels, NL **);
  nodelist->tails = DIRECTED ? Calloc(strat_nlevels, NL **) : nodelist->boths;
  nodelist->heads = DIRECTED ? Calloc(strat_nlevels, NL **) : nodelist->boths;
  
  for(int i = 0; i < strat_nlevels; i++) {
    nodelist->boths[i] = Calloc(bd_nlevels, NL *);
    if(DIRECTED) {
      nodelist->tails[i] = Calloc(bd_nlevels, NL *);
      nodelist->heads[i] = Calloc(bd_nlevels, NL *);
    }
    for(int j = 0; j < bd_nlevels; j++) {
      nodelist->boths[i][j] = NLInitialize(jointattrcounts[i*bd_nlevels + j], nodelist->bothpos);
      if(DIRECTED) {
        nodelist->tails[i][j] = NLInitialize(jointattrcounts[i*bd_nlevels + j], nodelist->tailpos);
        nodelist->heads[i][j] = NLInitialize(jointattrcounts[i*bd_nlevels + j], nodelist->headpos);
      }
    }
  }
  
  for(Vertex vertex = 1; vertex <= N_NODES; vertex++) {
    int strat_val = nodelist->strat_vattr[vertex];
    int bd_val = nodelist->bd_vattr[vertex];
    
    if(DIRECTED) {
      if(IN_DEG[vertex] < maxin && OUT_DEG[vertex] < maxout) {
        NLInsert(nodelist->boths[strat_val][bd_val], vertex);        
      } else if(OUT_DEG[vertex] < maxout) {
        NLInsert(nodelist->tails[strat_val][bd_val], vertex);
      } else if(IN_DEG[vertex] < maxin) {
        NLInsert(nodelist->heads[strat_val][bd_val], vertex);
      }        
    } else if(IN_DEG[vertex] + OUT_DEG[vertex] < maxout) {
      NLInsert(nodelist->boths[strat_val][bd_val], vertex);        
    }
  }

  return nodelist;
}

static inline Dyad NodeListDyadCount(NodeList *nodelist, int stratmixingtype) {
  NL **tailboths = nodelist->boths[nodelist->strat_tails[stratmixingtype]];
  NL **headboths = nodelist->boths[nodelist->strat_heads[stratmixingtype]];
  
  NL **tails = nodelist->tails[nodelist->strat_tails[stratmixingtype]];
  NL **heads = nodelist->heads[nodelist->strat_heads[stratmixingtype]];
  
  int *tailtypes = nodelist->bd_tails;
  int *headtypes = nodelist->bd_heads;

  int diagonal = nodelist->strat_diagonal[stratmixingtype];
  
  Dyad dyadcount = 0;
  for(int i = 0; i < nodelist->lengths[diagonal]; i++) {
    if(nodelist->directed) {
      dyadcount += NLDyadCount(tailboths[tailtypes[i]], headboths[headtypes[i]], diagonal && nodelist->bd_diagonal[i], TRUE);
      dyadcount += NLDyadCount(tailboths[tailtypes[i]], heads[headtypes[i]], FALSE, TRUE);
      dyadcount += NLDyadCount(tails[tailtypes[i]], headboths[headtypes[i]], FALSE, TRUE);
      dyadcount += NLDyadCount(tails[tailtypes[i]], heads[headtypes[i]], FALSE, TRUE);        
    } else {
      dyadcount += NLDyadCount(tails[tailtypes[i]], heads[headtypes[i]], diagonal && nodelist->bd_diagonal[i], nodelist->directed);
    }
  }
  return dyadcount;
}

static inline int NodeListDyadCountPositive(NodeList *nodelist, int stratmixingtype) {
  NL **tailboths = nodelist->boths[nodelist->strat_tails[stratmixingtype]];
  NL **headboths = nodelist->boths[nodelist->strat_heads[stratmixingtype]];
  
  NL **tails = nodelist->tails[nodelist->strat_tails[stratmixingtype]];
  NL **heads = nodelist->heads[nodelist->strat_heads[stratmixingtype]];
  
  int *tailtypes = nodelist->bd_tails;
  int *headtypes = nodelist->bd_heads;

  int diagonal = nodelist->strat_diagonal[stratmixingtype];
  
  for(int i = 0; i < nodelist->lengths[diagonal]; i++) {
    if(nodelist->directed) {
      if(NLDyadCountPositive(tailboths[tailtypes[i]], headboths[headtypes[i]], diagonal && nodelist->bd_diagonal[i]) ||
         NLDyadCountPositive(tailboths[tailtypes[i]], heads[headtypes[i]], FALSE) ||
         NLDyadCountPositive(tails[tailtypes[i]], headboths[headtypes[i]], FALSE) ||
         NLDyadCountPositive(tails[tailtypes[i]], heads[headtypes[i]], FALSE)) {
           return TRUE;
      }
    } else {
      if(NLDyadCountPositive(tails[tailtypes[i]], heads[headtypes[i]], diagonal && nodelist->bd_diagonal[i])) {
        return TRUE;
      }
    }
  }
  return FALSE;
}

static inline void NodeListToggleKnownIf(Vertex tail, Vertex head, NodeList *nodelist, int nodeflag, int tailcondition, int headcondition) {
  if(tailcondition) {
    if(nodelist->directed && (nodelist->bothpos[tail] || nodelist->headpos[tail])) {
      NLToggleKnown(nodelist->boths[nodelist->strat_vattr[tail]][nodelist->bd_vattr[tail]], tail, nodelist->bothpos[tail]);          
      NLToggleKnown(nodelist->heads[nodelist->strat_vattr[tail]][nodelist->bd_vattr[tail]], tail, nodelist->headpos[tail]);              
    } else {
      NLToggleKnown(nodelist->tails[nodelist->strat_vattr[tail]][nodelist->bd_vattr[tail]], tail, nodeflag);        
    }
  }

  if(headcondition) {
    if(nodelist->directed && (nodelist->bothpos[head] || nodelist->tailpos[head])) {
      NLToggleKnown(nodelist->boths[nodelist->strat_vattr[head]][nodelist->bd_vattr[head]], head, nodelist->bothpos[head]);          
      NLToggleKnown(nodelist->tails[nodelist->strat_vattr[head]][nodelist->bd_vattr[head]], head, nodelist->tailpos[head]);          
    } else {
      NLToggleKnown(nodelist->heads[nodelist->strat_vattr[head]][nodelist->bd_vattr[head]], head, nodeflag);
    }
  }
}

static inline void NodeListGetRandWithCount(Vertex *tail, Vertex *head, NodeList *nodelist, int stratmixingtype, Dyad dyadcount) {    
  NL **tailboths = nodelist->boths[nodelist->strat_tails[stratmixingtype]];
  NL **headboths = nodelist->boths[nodelist->strat_heads[stratmixingtype]];
  
  NL **tails = nodelist->tails[nodelist->strat_tails[stratmixingtype]];
  NL **heads = nodelist->heads[nodelist->strat_heads[stratmixingtype]];

  int *tailtypes = nodelist->bd_tails;
  int *headtypes = nodelist->bd_heads;

  int diagonal = nodelist->strat_diagonal[stratmixingtype];
                                              
  Dyad dyadindex = 2*dyadcount*unif_rand();

  for(int i = 0; i < nodelist->lengths[diagonal]; i++) {
    Dyad dyadsthistype;
    
    if(nodelist->directed) {
      dyadsthistype = NL2DyadCount(tailboths[tailtypes[i]], headboths[headtypes[i]], diagonal && nodelist->bd_diagonal[i], TRUE);
      if(dyadindex < dyadsthistype) {
        NLPut2Dyad(tail, head, tailboths[tailtypes[i]], headboths[headtypes[i]], dyadindex, diagonal && nodelist->bd_diagonal[i], TRUE);
        return;
      } else {
        dyadindex -= dyadsthistype;  
      }

      dyadsthistype = NL2DyadCount(tailboths[tailtypes[i]], heads[headtypes[i]], FALSE, TRUE);
      if(dyadindex < dyadsthistype) {
        NLPut2Dyad(tail, head, tailboths[tailtypes[i]], heads[headtypes[i]], dyadindex, FALSE, TRUE);
        return;
      } else {
        dyadindex -= dyadsthistype;  
      }

      dyadsthistype = NL2DyadCount(tails[tailtypes[i]], headboths[headtypes[i]], FALSE, TRUE);
      if(dyadindex < dyadsthistype) {
        NLPut2Dyad(tail, head, tails[tailtypes[i]], headboths[headtypes[i]], dyadindex, FALSE, TRUE);
        return;
      } else {
        dyadindex -= dyadsthistype;  
      }
      
      dyadsthistype = NL2DyadCount(tails[tailtypes[i]], heads[headtypes[i]], FALSE, TRUE);
      if(dyadindex < dyadsthistype) {
        NLPut2Dyad(tail, head, tails[tailtypes[i]], heads[headtypes[i]], dyadindex, FALSE, TRUE);
        return;
      } else {
        dyadindex -= dyadsthistype;  
      }
    } else {
      dyadsthistype = NL2DyadCount(tails[tailtypes[i]], heads[headtypes[i]], diagonal && nodelist->bd_diagonal[i], nodelist->directed);
      if(dyadindex < dyadsthistype) {
        NLPut2Dyad(tail, head, tails[tailtypes[i]], heads[headtypes[i]], dyadindex, diagonal && nodelist->bd_diagonal[i], nodelist->directed);
        return;
      } else {
        dyadindex -= dyadsthistype;  
      }        
    }
  }
}

static inline void NodeListGetRand(Vertex *tail, Vertex *head, NodeList *nodelist, int stratmixingtype) {
  NodeListGetRandWithCount(tail, head, nodelist, stratmixingtype, NodeListDyadCount(nodelist, stratmixingtype));
}

static inline void NodeListDestroy(NodeList *nodelist) {
  for(int i = 0; i < nodelist->strat_nlevels; i++) {
    for(int j = 0; j < nodelist->bd_nlevels; j++) {
      NLDestroy(nodelist->boths[i][j]);
      if(nodelist->directed) {
        NLDestroy(nodelist->tails[i][j]);
        NLDestroy(nodelist->heads[i][j]);
      }
    }
  }
  
  Free(nodelist->bothpos);
  if(nodelist->directed) {
    Free(nodelist->tailpos);
    Free(nodelist->headpos);
  }
  
  Free(nodelist->strat_diagonal);
  Free(nodelist->bd_diagonal);
  
  Free(nodelist);
}

static inline Dyad NodeListDyadCountOnToggle(Vertex tail, Vertex head, NodeList *nodelist, int stratmixingtype, int change, int tailcondition, int headcondition) {
  if(tailcondition) {
    if(nodelist->directed && (nodelist->bothpos[tail] || nodelist->headpos[tail])) {
      nodelist->boths[nodelist->strat_vattr[tail]][nodelist->bd_vattr[tail]]->length += change;          
      nodelist->heads[nodelist->strat_vattr[tail]][nodelist->bd_vattr[tail]]->length -= change;              
    } else {
      nodelist->tails[nodelist->strat_vattr[tail]][nodelist->bd_vattr[tail]]->length += change;        
    }
  }

  if(headcondition) {
    if(nodelist->directed && (nodelist->bothpos[head] || nodelist->tailpos[head])) {
      nodelist->boths[nodelist->strat_vattr[head]][nodelist->bd_vattr[head]]->length += change;        
      nodelist->tails[nodelist->strat_vattr[head]][nodelist->bd_vattr[head]]->length -= change;        
    } else {
      nodelist->heads[nodelist->strat_vattr[head]][nodelist->bd_vattr[head]]->length += change;
    }
  }

  Dyad dyadcount = NodeListDyadCount(nodelist, stratmixingtype);

  if(tailcondition) {
    if(nodelist->directed && (nodelist->bothpos[tail] || nodelist->headpos[tail])) {
      nodelist->boths[nodelist->strat_vattr[tail]][nodelist->bd_vattr[tail]]->length -= change;          
      nodelist->heads[nodelist->strat_vattr[tail]][nodelist->bd_vattr[tail]]->length += change;              
    } else {
      nodelist->tails[nodelist->strat_vattr[tail]][nodelist->bd_vattr[tail]]->length -= change;        
    }
  }

  if(headcondition) {
    if(nodelist->directed && (nodelist->bothpos[head] || nodelist->tailpos[head])) {
      nodelist->boths[nodelist->strat_vattr[head]][nodelist->bd_vattr[head]]->length -= change;        
      nodelist->tails[nodelist->strat_vattr[head]][nodelist->bd_vattr[head]]->length += change;        
    } else {
      nodelist->heads[nodelist->strat_vattr[head]][nodelist->bd_vattr[head]]->length -= change;
    }
  }

  return dyadcount;
}


#endif
