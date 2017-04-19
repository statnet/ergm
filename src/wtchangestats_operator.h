#ifndef WTCHANGESTATS_OPERATOR_H
#define WTCHANGESTATS_OPERATOR_H

#define STRICT_Wt_HEADERS
#include "changestat_operator.h"
#include "wtchangestat_operator.h"
#include "ergm_model.h"
#include "ergm_wtchangestat.h"
#include "ergm_wtmodel.h"
#include "ergm_storage.h"

typedef struct{Network nw; Model *m;} StoreNetAndModel;
typedef struct{Network nw; WtModel *m;} StoreNetAndWtModel;

#endif // WTCHANGESTATS_OPERATOR_H
