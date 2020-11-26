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

#endif
