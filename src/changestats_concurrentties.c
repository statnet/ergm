#include "changestats_concurrentties.h"

/*****************
 changestat: d_concurrent_ties
*****************/
D_CHANGESTAT_FN(d_concurrent_ties) { 
  int i, echange;
  Vertex h, t, hdeg, tdeg;

  CHANGE_STAT[0] = 0.0;  
  FOR_EACH_TOGGLE(i) {
    h = heads[i];
    t = tails[i];
    echange = IS_OUTEDGE(h, t) ? -1 : 1;
    hdeg = OUT_DEG[h];
    tdeg = IN_DEG[t];
    if(!DIRECTED){
      hdeg += IN_DEG[h];
      tdeg += OUT_DEG[t];
    }
    if(echange>0){
      // When adding ties, only count those beyond the node's first.
      if(hdeg>=1) CHANGE_STAT[0]++;
      if(tdeg>=1) CHANGE_STAT[0]++;
    }else{
      // When removing ties, only count if the second tie or beyond is removed.
      if(hdeg>=2) CHANGE_STAT[0]--;
      if(tdeg>=2) CHANGE_STAT[0]--;
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
  int i, j, echange, hattr, tattr;
  Vertex h, t, hdeg, tdeg;

  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    h = heads[i];
    t = tails[i];
    echange = IS_OUTEDGE(h, t) ? -1 : 1;
    hdeg = OUT_DEG[h];
    tdeg = IN_DEG[t];
    if(!DIRECTED){
      hdeg += IN_DEG[h];
      tdeg += OUT_DEG[t];
    }
    hattr = INPUT_PARAM[N_CHANGE_STATS + h - 1]; 
    tattr = INPUT_PARAM[N_CHANGE_STATS + t - 1];
    for(j = 0; j < N_CHANGE_STATS; j++) {
      if(echange>0){
	// When adding ties, only count those beyond the node's first.
	if(hattr == INPUT_PARAM[j] && hdeg>=1) CHANGE_STAT[j]++;
	if(tattr == INPUT_PARAM[j] && tdeg>=1) CHANGE_STAT[j]++;
      }else{
	// When removing ties, only count if the second tie or beyond is removed.
	if(hattr == INPUT_PARAM[j] && hdeg>=2) CHANGE_STAT[j]--;
	if(tattr == INPUT_PARAM[j] && tdeg>=2) CHANGE_STAT[j]--;
      }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}


