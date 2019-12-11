#ifndef _ERGM_WTSTATE_H_
#define _ERGM_WTSTATE_H_

#include "ergm_wtedgetree.h"
#include "ergm_wtchangestat.h"
#include "ergm_wtMHproposal.h"
#include "ergm_wtmodel.h"

#define ARGS_WTSTATE SEXP stateR
#define ARGS_WTLASTTOGGLE SEXP time, SEXP lasttoggle

#define YES_WTSTATE stateR, FALSE, FALSE
#define YES_WTSTATE_EMPTY_NO_INIT_S stateR, TRUE, TRUE
#define YES_WTLASTTOGGLE length(time)>0, asInteger(time), INTEGER(lasttoggle)

#define NO_WTLASTTOGGLE FALSE, 0, NULL

typedef struct{
  double *stats;
  WtNetwork *nwp;
  WtModel *m;
  WtMHProposal *MHp;
} WtErgmState;

WtErgmState *WtErgmStateInit(SEXP stateR,
                             Rboolean empty, Rboolean noinit_s,
                             // Network state
                             Rboolean timings, int time, int *lasttoggle);
SEXP WtErgmStateRSave(SEXP startR, WtErgmState *s);
void WtErgmStateDestroy(WtErgmState *s);

#endif
