#ifndef _ERGM_WTSTATE_H_
#define _ERGM_WTSTATE_H_

#include "ergm_wtedgetree.h"
#include "ergm_wtchangestat.h"
#include "ergm_wtMHproposal.h"
#include "ergm_wtmodel.h"

#define ARGS_WTSTATE SEXP stateR

#define YES_WTSTATE stateR, FALSE, FALSE
#define YES_WTSTATE_EMPTY_NO_INIT_S stateR, TRUE, TRUE

typedef struct{
  double *stats;
  WtNetwork *nwp;
  WtModel *m;
  WtMHProposal *MHp;
} WtErgmState;

WtErgmState *WtErgmStateInit(SEXP stateR,
                             Rboolean empty, Rboolean noinit_s);
SEXP WtErgmStateRSave(SEXP startR, WtErgmState *s);
void WtErgmStateDestroy(WtErgmState *s);

#endif
