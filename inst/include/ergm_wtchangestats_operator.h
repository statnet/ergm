/*  File inst/include/ergm_wtchangestats_operator.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2023 Statnet Commons
 */
#ifndef _ERGM_WTCHANGESTATS_OPERATOR_H_
#define _ERGM_WTCHANGESTATS_OPERATOR_H_

#define STRICT_Wt_HEADERS
#include "ergm_wtmodel.h"

typedef struct{WtModel *m; double *stats;} StoreWtModelAndStats;

#endif // _ERGM_WTCHANGESTATS_OPERATOR_H_
