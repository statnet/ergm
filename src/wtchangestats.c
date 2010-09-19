#include "wtchangestats.h"

/********************  changestats:   N    ***********/

/*****************
 changestat: d_nonzero
*****************/
WtD_CHANGESTAT_FN(d_nonzero) {
  double oldwt;
  int i;
  
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i){
    GETOLDWT(i);
    CHANGE_STAT[0] += (weights[i]!=0) - (OLDWT!=0);
    SETWT_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_SETWTS(i);
}



/********************  changestats:   S    ***********/

/*****************
 changestat: d_sum
*****************/
WtD_CHANGESTAT_FN(d_sum) {
  double oldwt;
  int i;
  
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i){
    GETOLDWT(i);
    CHANGE_STAT[0] += weights[i]-OLDWT;
    SETWT_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_SETWTS(i);
}


