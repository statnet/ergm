#include "wtchangestat.h"

D_CHANGESTAT_FN(d_from_s) { 
  int i;
  double current, OLDWT;

  (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */
  current = mtp->dstats[0];
  FOR_EACH_TOGGLE(i) { SETWT_WITH_BACKUP(i) }
  (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */
  mtp->dstats[0] -= current;
  FOR_EACH_TOGGLE(i) { UNDO_SETWT(i) }
}
