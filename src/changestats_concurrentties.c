/*  File src/changestats_concurrentties.c in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2021 Statnet Commons
 */
#include "ergm_edgetree.h"
#include "ergm_changestat.h"

/*****************
 changestat: d_concurrent_ties
*****************/
C_CHANGESTAT_FN(c_concurrent_ties) { 
  int echange;
  Vertex taildeg, headdeg;

    echange = IS_OUTEDGE(tail, head) ? -1 : 1;
    taildeg = OUT_DEG[tail];
    headdeg = IN_DEG[head];
    if(!DIRECTED){
      taildeg += IN_DEG[tail];
      headdeg += OUT_DEG[head];
    }
    if(echange>0){
      // When adding ties, only count those beyond the node's first.
      if(taildeg>=1) CHANGE_STAT[0]++;
      if(headdeg>=1) CHANGE_STAT[0]++;
    }else{
      // When removing ties, only count if the second tie or beyond is removed.
      if(taildeg>=2) CHANGE_STAT[0]--;
      if(headdeg>=2) CHANGE_STAT[0]--;
    }
}

/*****************
 changestat: d_concurrent_ties_by_attr
*****************/
C_CHANGESTAT_FN(c_concurrent_ties_by_attr) { 
  /*The inputparams are assumed to be set up as follows:
    The first 2*nstats values are in pairs:  (degree, attrvalue)
    The values following the first 2*nstats values are the nodal attributes.
  */
  int j, echange, tailattr, headattr;
  Vertex taildeg, headdeg;

    echange = IS_OUTEDGE(tail, head) ? -1 : 1;
    taildeg = OUT_DEG[tail];
    headdeg = IN_DEG[head];
    if(!DIRECTED){
      taildeg += IN_DEG[tail];
      headdeg += OUT_DEG[head];
    }
    tailattr = INPUT_PARAM[N_CHANGE_STATS + tail - 1]; 
    headattr = INPUT_PARAM[N_CHANGE_STATS + head - 1];
    for(j = 0; j < N_CHANGE_STATS; j++) {
      if(echange>0){
        // When adding ties, only count those beyond the node's first.
        if(tailattr == INPUT_PARAM[j] && taildeg>=1) CHANGE_STAT[j]++;
        if(headattr == INPUT_PARAM[j] && headdeg>=1) CHANGE_STAT[j]++;
      }else{
        // When removing ties, only count if the second tie or beyond is removed.
        if(tailattr == INPUT_PARAM[j] && taildeg>=2) CHANGE_STAT[j]--;
        if(headattr == INPUT_PARAM[j] && headdeg>=2) CHANGE_STAT[j]--;
      }
    }
}


