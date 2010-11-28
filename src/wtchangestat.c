#include "wtchangestat.h"

WtD_CHANGESTAT_FN(d_from_s) { 
  double current;

  (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */
  current = mtp->dstats[0];
  // Note: This cannot be abstracted into EXEC_THROUGH_TOGGLES.
  FOR_EACH_TOGGLE() { GETTOGGLEINFO(); SETWT_WITH_BACKUP(); }
  (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */
  mtp->dstats[0] -= current;
  FOR_EACH_TOGGLE() { UNDO_SETWT(); }
}
