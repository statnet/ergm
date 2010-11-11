#ifndef WTCHANGESTATS_H
#define WTCHANGESTATS_H

#include "wtedgetree.h"
#include "wtchangestat.h"

/********************  changestats:   C    ***********/
WtD_CHANGESTAT_FN(d_cyclicweighs_max); WtS_CHANGESTAT_FN(s_cyclicweights_max);

/********************  changestats:   N    ***********/
WtD_CHANGESTAT_FN(d_nonzero);

/********************  changestats:   S    ***********/
WtD_CHANGESTAT_FN(d_sum);

/********************  changestats:   T    ***********/
WtD_CHANGESTAT_FN(d_transitiveweighs_max); WtS_CHANGESTAT_FN(s_transitiveweights_max);

              
#endif
