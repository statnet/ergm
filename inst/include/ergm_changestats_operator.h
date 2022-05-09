/*  File inst/include/ergm_changestats_operator.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2022 Statnet Commons
 */
#ifndef _ERGM_CHANGESTATS_OPERATOR_H_
#define _ERGM_CHANGESTATS_OPERATOR_H_

#include "ergm_model.h"

typedef struct{Model *m; double *stats;} StoreModelAndStats;

#endif // _ERGM_CHANGESTATS_OPERATOR_H_
