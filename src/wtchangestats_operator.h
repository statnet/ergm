#ifndef _WTCHANGESTATS_OPERATOR_H_
#define _WTCHANGESTATS_OPERATOR_H_

#define STRICT_Wt_HEADERS
#include "ergm_changestat_operator.h"
#include "ergm_wtchangestat_operator.h"
#include "ergm_model.h"
#include "ergm_wtchangestat.h"
#include "ergm_wtmodel.h"
#include "ergm_storage.h"

typedef struct{Network nw; Model *m;} StoreNetAndModel;
typedef struct{Network nw; WtModel *m;} StoreNetAndWtModel;

#endif // WTCHANGESTATS_OPERATOR_H
