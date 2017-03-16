#include "changestats_indices.h"

/*****************                       
 changestat: d_indices
*****************/
C_CHANGESTAT_FN(c_indices) { 
  ZERO_ALL_CHANGESTATS(i);
    int edgeflag = IS_OUTEDGE(tail,head);
    CHANGE_STAT[0] += edgeflag ? -tail : tail;
    CHANGE_STAT[1] += edgeflag ? -head : head;
}
