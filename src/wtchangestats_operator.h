/*  File src/wtchangestats_operator.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2023 Statnet Commons
 */
#ifndef _WTCHANGESTATS_OPERATOR_H_
#define _WTCHANGESTATS_OPERATOR_H_

#define STRICT_Wt_HEADERS
#include "ergm_wtchangestat_operator.h"
#include "ergm_changestat_operator.h"
#include "ergm_model.h"
#include "ergm_wtmodel.h"
#include "ergm_model.h"
#include "ergm_wtchangestat.h"
#include "ergm_wtmodel.h"
#include "ergm_storage.h"

typedef struct{Network *nwp; Model *m;} StoreNetAndModel;
typedef struct{Network *nwp; WtModel *m;} StoreNetAndWtModel;

#endif // WTCHANGESTATS_OPERATOR_H
