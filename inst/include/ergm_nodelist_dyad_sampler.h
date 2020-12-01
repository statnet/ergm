#ifndef _ERGM_NODELIST_DYAD_SAMPLER_H_
#define _ERGM_NODELIST_DYAD_SAMPLER_H_

#include "ergm_edgetree_types.h"

// a node-list dyad sampler that is general enough
// to cover all the new BD/Strat proposals

static inline void GetRandDyadFromLists(Vertex *tail, 
                                        Vertex *head, 
                                        Vertex **tails,
                                        Vertex **heads,
                                        int *tailattrs,
                                        int *headattrs,
                                        int *tailcounts, 
                                        int *headcounts, 
                                        int length, 
                                        Dyad dyadcount,
                                        int diagonal, 
                                        int directed) {
                                            
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

// for simplicity of U_FNs, adding this trivial updating function as well
static inline void NodeListToggleKnown(Vertex node, Vertex *nodevec, int *nodepos, int *length, int nodeflag) {
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

static inline Dyad NodeListDyadCount(int *tailcounts, int *headcounts, int *tailtypes, int *headtypes, int length, int diagonal, int directed) {
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

static inline int NodeListDyadCountPositive(int *tailcounts, int *headcounts, int *tailtypes, int *headtypes, int length, int diagonal) {
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

#endif
