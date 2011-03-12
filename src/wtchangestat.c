#include "wtchangestat.h"

WtD_CHANGESTAT_FN(d_from_s) { 
  (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */
  memcpy(mtp->statcache,mtp->dstats,N_CHANGE_STATS*sizeof(double));
  // Note: This cannot be abstracted into EXEC_THROUGH_TOGGLES.
  FOR_EACH_TOGGLE() { GETTOGGLEINFO(); SETWT_WITH_BACKUP(); }
  (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */
  for(unsigned int i=0; i<N_CHANGE_STATS; i++) mtp->dstats[i] -= mtp->statcache[i];
  FOR_EACH_TOGGLE() { UNDO_SETWT(); }
}
