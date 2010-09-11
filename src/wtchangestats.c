#include "wtchangestats.h"

/********************  changestats:   S    ***********/

/*****************
 changestat: d_sum
*****************/
WtD_CHANGESTAT_FN(d_sum) {
  double oldwt;
  int i;
  
  CHANGE_STAT[0] = 0.0;
  for (i=0; i < ntoggles; i++)
    {
      GETOLDWT(i);
      CHANGE_STAT[0] += weights[i]-OLDWT;
      SETWT_IF_MORE_TO_COME(i);
    }
  UNDO_PREVIOUS_SETWTS(i);
}


