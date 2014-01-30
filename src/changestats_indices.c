#include "changestats_indices.h"

/*****************                       
 changestat: d_indices
*****************/
D_CHANGESTAT_FN(d_indices) { 
  int i;

  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    Vertex t = tails[i];
    Vertex h = heads[i]; 
    int edgeflag = IS_OUTEDGE(t,h);
    CHANGE_STAT[0] += edgeflag ? -t : t;
    CHANGE_STAT[1] += edgeflag ? -h : h;
    TOGGLE_IF_MORE_TO_COME(i); /* Needed in case of multiple toggles */
  }
  UNDO_PREVIOUS_TOGGLES(i); /* Needed on exit in case of multiple toggles */
}
