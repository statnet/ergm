#ifndef _ERGM_NODELIST_H_
#define _ERGM_NODELIST_H_

#include "ergm_MHproposal.h"
#include "ergm_edgetree_types.h"
#include "ergm_changestat.h"

typedef struct {
  Vertex ***nodes;
  int **attrcounts;
  int *nodepos;
  
  int *strat_tails;
  int *strat_heads;
  
  int *bd_tails;
  int *bd_heads;

  int *strat_vattr;
  int *bd_vattr;
  
  int *lengths;
  
  int directed;

  int strat_nlevels;
  int bd_nlevels;
} NodeList;

static inline NodeList *NodeListInitialize(int bound, int *strat_vattr, int strat_nlevels, int strat_nmixtypes, int *strat_tails, int *strat_heads, int *strat_attrcounts,
                                           int *bd_vattr, int bd_nlevels, int nondiag_bd_nmixtypes, int diag_bd_nmixtypes, int *bd_tails, int *bd_heads, int *bd_attrcounts, 
                                           int *jointattrcounts, Network *nwp) {
  NodeList *nodelist = Calloc(1, NodeList);
  nodelist->directed = DIRECTED;
  
  if(!strat_vattr) {
    strat_nlevels = 1;
    nodelist->strat_tails = Calloc(1, int);
    nodelist->strat_heads = Calloc(1, int);
    nodelist->strat_vattr = Calloc(N_NODES, int);
  } else {
    nodelist->strat_tails = Calloc(strat_nmixtypes, int);
    nodelist->strat_heads = Calloc(strat_nmixtypes, int);
    nodelist->strat_vattr = Calloc(N_NODES, int);
    memcpy(nodelist->strat_tails, strat_tails, strat_nmixtypes*sizeof(int));
    memcpy(nodelist->strat_heads, strat_heads, strat_nmixtypes*sizeof(int));
    memcpy(nodelist->strat_vattr, strat_vattr, N_NODES*sizeof(int));
  }
  nodelist->strat_vattr--; // decrement so node indices line up correctly
  
  nodelist->lengths = Calloc(2, int);
  if(!bd_vattr) {
    bd_nlevels = 1;
    nodelist->bd_tails = Calloc(1, int);
    nodelist->bd_heads = Calloc(1, int);
    nodelist->bd_vattr = Calloc(N_NODES, int);
    nodelist->lengths[0] = 1;
    nodelist->lengths[1] = 1;
  } else {
    nodelist->bd_tails = Calloc(nondiag_bd_nmixtypes, int);
    nodelist->bd_heads = Calloc(nondiag_bd_nmixtypes, int);
    nodelist->bd_vattr = Calloc(N_NODES, int);
    memcpy(nodelist->bd_tails, bd_tails, nondiag_bd_nmixtypes*sizeof(int));
    memcpy(nodelist->bd_heads, bd_heads, nondiag_bd_nmixtypes*sizeof(int));
    memcpy(nodelist->bd_vattr, bd_vattr, N_NODES*sizeof(int));
    nodelist->lengths[0] = nondiag_bd_nmixtypes;
    nodelist->lengths[1] = diag_bd_nmixtypes;
  }
  nodelist->bd_vattr--; // decrement so node indices line up correctly
  
  Vertex ***nodes = Calloc(strat_nlevels, Vertex **);
  for(int i = 0; i < strat_nlevels; i++) {
    nodes[i] = Calloc(bd_nlevels, Vertex *);
    for(int j = 0; j < bd_nlevels; j++) {
      nodes[i][j] = Calloc(jointattrcounts ? jointattrcounts[i*bd_nlevels + j] : strat_attrcounts ? strat_attrcounts[i] : bd_attrcounts[j], Vertex);
    }
  }
  
  int **attrcounts = Calloc(strat_nlevels, int *);
  for(int i = 0; i < strat_nlevels; i++) {
    attrcounts[i] = Calloc(bd_nlevels, int);
  }
  
  int *nodepos = Calloc(N_NODES, int);
  nodepos--;
  
  for(Vertex vertex = 1; vertex <= N_NODES; vertex++) {
    if(IN_DEG[vertex] + OUT_DEG[vertex] < bound) {
      int strat_val = nodelist->strat_vattr[vertex];
      int bd_val = nodelist->bd_vattr[vertex];
      nodes[strat_val][bd_val][attrcounts[strat_val][bd_val]] = vertex;
      nodepos[vertex] = attrcounts[strat_val][bd_val];
      attrcounts[strat_val][bd_val]++;
    }
  }
  
  nodelist->nodes = nodes;
  nodelist->attrcounts = attrcounts;
  nodelist->nodepos = nodepos;
  
  nodelist->strat_nlevels = strat_nlevels;
  nodelist->bd_nlevels = bd_nlevels;
  
  return nodelist;
}

static inline Dyad NodeListDyadCount(NodeList *nodelist, int stratmixingtype) {
  int *tailcounts = nodelist->attrcounts[nodelist->strat_tails[stratmixingtype]];
  int *headcounts = nodelist->attrcounts[nodelist->strat_heads[stratmixingtype]];
  int *tailtypes = nodelist->bd_tails;
  int *headtypes = nodelist->bd_heads;
  int diagonal = nodelist->strat_tails[stratmixingtype] == nodelist->strat_heads[stratmixingtype];
  int length = nodelist->lengths[diagonal];
  int directed = nodelist->directed;

  Dyad dyadcount = 0;
  for(int i = 0; i < length; i++) {
    if(diagonal && tailtypes[i] == headtypes[i]) {
      if(directed) {
        dyadcount += (Dyad)tailcounts[tailtypes[i]]*(headcounts[headtypes[i]] - 1);
      } else {
        dyadcount += (Dyad)tailcounts[tailtypes[i]]*(headcounts[headtypes[i]] - 1)/2;
      }
    } else {
      dyadcount += (Dyad)tailcounts[tailtypes[i]]*headcounts[headtypes[i]];
    }
  }
  return dyadcount;
}

static inline int NodeListDyadCountPositive(NodeList *nodelist, int stratmixingtype) {
  int *tailcounts = nodelist->attrcounts[nodelist->strat_tails[stratmixingtype]];
  int *headcounts = nodelist->attrcounts[nodelist->strat_heads[stratmixingtype]];
  int *tailtypes = nodelist->bd_tails;
  int *headtypes = nodelist->bd_heads;
  int diagonal = nodelist->strat_tails[stratmixingtype] == nodelist->strat_heads[stratmixingtype];
  int length = nodelist->lengths[diagonal];

  for(int i = 0; i < length; i++) {
    if(diagonal && tailtypes[i] == headtypes[i]) {
      if(tailcounts[tailtypes[i]] > 1) {
        return TRUE;
      }
    } else {
      if(tailcounts[tailtypes[i]] > 0 && headcounts[headtypes[i]] > 0) {
        return TRUE;
      }
    }
  }
  return FALSE;
}

static inline void NodeListToggleKnown(Vertex node, NodeList *nodelist, int nodeflag) {
  Vertex *nodevec = nodelist->nodes[nodelist->strat_vattr[node]][nodelist->bd_vattr[node]];
  int *nodepos = nodelist->nodepos;
  int *length = &nodelist->attrcounts[nodelist->strat_vattr[node]][nodelist->bd_vattr[node]];
  
  if(nodeflag) { // present; remove node from the list
    nodevec[nodepos[node]] = nodevec[*length - 1];  // move last node into node's position
    nodepos[nodevec[*length - 1]] = nodepos[node];  // update formerly last node's position
    (*length)--;  // update length
  } else { // absent; add node to the list
    nodevec[*length] = node;  // place node in last position
    nodepos[node] = *length;  // set position for node
    (*length)++;  // update length  
  }
}

static inline void NodeListToggleKnownIf(Vertex tail, Vertex head, NodeList *nodelist, int nodeflag, int tailcondition, int headcondition) {
  if(tailcondition) {
    NodeListToggleKnown(tail, nodelist, nodeflag);
  }
  if(headcondition) {
    NodeListToggleKnown(head, nodelist, nodeflag);
  }
}

static inline void NodeListGetRandWithCount(Vertex *tail, Vertex *head, NodeList *nodelist, int stratmixingtype, Dyad dyadcount) {    
  Vertex **tails = nodelist->nodes[nodelist->strat_tails[stratmixingtype]];
  Vertex **heads = nodelist->nodes[nodelist->strat_heads[stratmixingtype]];

  int *tailattrs = nodelist->bd_tails;
  int *headattrs = nodelist->bd_heads;

  int *tailcounts = nodelist->attrcounts[nodelist->strat_tails[stratmixingtype]];
  int *headcounts = nodelist->attrcounts[nodelist->strat_heads[stratmixingtype]];

  int diagonal = nodelist->strat_tails[stratmixingtype] == nodelist->strat_heads[stratmixingtype];
  int length = nodelist->lengths[diagonal];

  int directed = nodelist->directed;
                                              
  Dyad dyadindex = 2*dyadcount*unif_rand();

  for(int i = 0; i < length; i++) {
    Dyad dyadsthistype;

    int tailcount = tailcounts[tailattrs[i]];
    int headcount = headcounts[headattrs[i]];

    // doubled for efficiency below
    if(diagonal && tailattrs[i] == headattrs[i]) {
      if(directed) {
        dyadsthistype = (Dyad)2*tailcount*(headcount - 1);                    
      } else {
        dyadsthistype = (Dyad)tailcount*(headcount - 1);          
      }
    } else {
      dyadsthistype = (Dyad)2*tailcount*headcount;
    }
    
    if(dyadindex < dyadsthistype) {
      int tailindex;
      int headindex;
      
      if(diagonal && tailattrs[i] == headattrs[i]) {
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

      if(tails[tailattrs[i]][tailindex] < heads[headattrs[i]][headindex] || directed) {
        *tail = tails[tailattrs[i]][tailindex];
        *head = heads[headattrs[i]][headindex];
      } else {
        *tail = heads[headattrs[i]][headindex];
        *head = tails[tailattrs[i]][tailindex];          
      }
      
      return;
    } else {
      dyadindex -= dyadsthistype;
    }
  } 
}

static inline void NodeListGetRand(Vertex *tail, Vertex *head, NodeList *nodelist, int stratmixingtype) {
  NodeListGetRandWithCount(tail, head, nodelist, stratmixingtype, NodeListDyadCount(nodelist, stratmixingtype));
}

static inline void NodeListDestroy(NodeList *nodelist) {
  for(int i = 0; i < nodelist->strat_nlevels; i++) {
    for(int j = 0; j < nodelist->bd_nlevels; j++) {
      Free(nodelist->nodes[i][j]);
    }
    Free(nodelist->nodes[i]);
  }
  Free(nodelist->nodes);
  
  for(int i = 0; i < nodelist->strat_nlevels; i++) {
    Free(nodelist->attrcounts[i]);
  }
  Free(nodelist->attrcounts);

  nodelist->nodepos++;
  Free(nodelist->nodepos);  
  Free(nodelist->strat_tails);  
  Free(nodelist->strat_heads);  
  Free(nodelist->bd_tails);  
  Free(nodelist->bd_heads);  
  nodelist->strat_vattr++;
  Free(nodelist->strat_vattr);
  nodelist->bd_vattr++;
  Free(nodelist->bd_vattr);
  Free(nodelist->lengths);  
  Free(nodelist);
}


static inline Dyad NodeListDyadCountOnToggle(Vertex tail, Vertex head, NodeList *nodelist, int stratmixingtype, int change, int tailcondition, int headcondition) {
  int *tailcounts = nodelist->attrcounts[nodelist->strat_tails[stratmixingtype]];
  int *headcounts = nodelist->attrcounts[nodelist->strat_heads[stratmixingtype]];
  int *tailtypes = nodelist->bd_tails;
  int *headtypes = nodelist->bd_heads;
  int diagonal = nodelist->strat_tails[stratmixingtype] == nodelist->strat_heads[stratmixingtype];
  int length = nodelist->lengths[diagonal];
  int directed = nodelist->directed;
  
  if(tailcondition) {
    nodelist->attrcounts[nodelist->strat_vattr[tail]][nodelist->bd_vattr[tail]] += change;
  }
  if(headcondition) {
    nodelist->attrcounts[nodelist->strat_vattr[head]][nodelist->bd_vattr[head]] += change;
  }
  Dyad dyadcount = 0;
  for(int i = 0; i < length; i++) {
    if(diagonal && tailtypes[i] == headtypes[i]) {
      if(directed) {
        dyadcount += (Dyad)tailcounts[tailtypes[i]]*(headcounts[headtypes[i]] - 1);
      } else {
        dyadcount += (Dyad)tailcounts[tailtypes[i]]*(headcounts[headtypes[i]] - 1)/2;
      }
    } else {
      dyadcount += (Dyad)tailcounts[tailtypes[i]]*headcounts[headtypes[i]];
    }
  }
  if(tailcondition) {
    nodelist->attrcounts[nodelist->strat_vattr[tail]][nodelist->bd_vattr[tail]] -= change;
  }
  if(headcondition) {
    nodelist->attrcounts[nodelist->strat_vattr[head]][nodelist->bd_vattr[head]] -= change;
  }
  return dyadcount;
}


#endif
