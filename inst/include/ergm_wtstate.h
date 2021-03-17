#ifndef _ERGM_WTSTATE_H_
#define _ERGM_WTSTATE_H_

#include "ergm_constants.h"
#include "ergm_wtedgetree.h"
#include "ergm_wtchangestat.h"
#include "ergm_wtMHproposal.h"
#include "ergm_wtmodel.h"

typedef struct{
  double *stats;
  WtNetwork *nwp;
  WtModel *m;
  WtMHProposal *MHp;
} WtErgmState;

WtErgmState *WtErgmStateInit(SEXP stateR,
                             unsigned int flags);
SEXP WtErgmStateRSave(SEXP startR, WtErgmState *s);
void WtErgmStateDestroy(WtErgmState *s);

#endif
