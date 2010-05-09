#include "changestats.h"

/*****************
 changestat: d_degcrossprod
*****************/
D_CHANGESTAT_FN(d_degcrossprod) { 
  int i, echange;
  Vertex h, t, hdeg, tdeg, node3;
  Edge e;

  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i) {
    h = heads[i];
    t = tails[i];
    echange = IS_OUTEDGE(h, t) ? -1 : 1;
    hdeg = OUT_DEG[h] + IN_DEG[h];
    tdeg = OUT_DEG[t] + IN_DEG[t];
    if(echange==1){
     CHANGE_STAT[0] += (hdeg + 1)*(tdeg + 1);
     STEP_THROUGH_OUTEDGES(t, e, node3) { /* step through outedges of tail */
      if(node3!=h) CHANGE_STAT[0] += (OUT_DEG[node3] + IN_DEG[node3]);
     }
     STEP_THROUGH_INEDGES(t, e, node3) { /* step through inedges of tail */
      if(node3!=h) CHANGE_STAT[0] += (OUT_DEG[node3] + IN_DEG[node3]);
     }
     STEP_THROUGH_OUTEDGES(h, e, node3) { /* step through outedges of head */
      if(node3!=t) CHANGE_STAT[0] += (OUT_DEG[node3] + IN_DEG[node3]);
     }
     STEP_THROUGH_INEDGES(h, e, node3) { /* step through inedges of head */
      if(node3!=t) CHANGE_STAT[0] += (OUT_DEG[node3] + IN_DEG[node3]);
     }
    }else{
     CHANGE_STAT[0] -= (hdeg)*(tdeg);
     STEP_THROUGH_OUTEDGES(t, e, node3) { /* step through outedges of tail */
      if(node3!=h) CHANGE_STAT[0] -= (OUT_DEG[node3] + IN_DEG[node3]);
     }
     STEP_THROUGH_INEDGES(t, e, node3) { /* step through inedges of tail */
      if(node3!=h) CHANGE_STAT[0] -= (OUT_DEG[node3] + IN_DEG[node3]);
     }
     STEP_THROUGH_OUTEDGES(h, e, node3) { /* step through outedges of head */
      if(node3!=t) CHANGE_STAT[0] -= (OUT_DEG[node3] + IN_DEG[node3]);
     }
     STEP_THROUGH_INEDGES(h, e, node3) { /* step through inedges of head */
      if(node3!=t) CHANGE_STAT[0] -= (OUT_DEG[node3] + IN_DEG[node3]);
     }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}
