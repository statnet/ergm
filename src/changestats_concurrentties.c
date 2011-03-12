#include "changestats_concurrentties.h"

/*****************
 changestat: d_concurrent_ties
*****************/
D_CHANGESTAT_FN(d_concurrent_ties) { 
  int i, echange;
  Vertex tail, head, taildeg, headdeg;

  CHANGE_STAT[0] = 0.0;  
  FOR_EACH_TOGGLE(i) {
    tail = TAIL(i);
    head = HEAD(i);
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
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_concurrent_ties_by_attr
*****************/
D_CHANGESTAT_FN(d_concurrent_ties_by_attr) { 
  /*The inputparams are assumed to be set up as follows:
    The first 2*nstats values are in pairs:  (degree, attrvalue)
    The values following the first 2*nstats values are the nodal attributes.
  */
  int i, j, echange, tailattr, headattr;
  Vertex tail, head, taildeg, headdeg;

  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    tail = TAIL(i);
    head = HEAD(i);
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
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}


