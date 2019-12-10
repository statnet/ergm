#ifndef _ERGM_WTSTATE_H_
#define _ERGM_WTSTATE_H_

#include "ergm_wtedgetree.h"
#include "ergm_wtchangestat.h"
#include "ergm_wtMHproposal.h"
#include "ergm_wtmodel.h"

#define ARGS_WTNW SEXP stateR
#define ARGS_WTMODEL SEXP mR
#define ARGS_WTMHPROPOSAL SEXP pR
#define ARGS_WTLASTTOGGLE SEXP time, SEXP lasttoggle

#define YES_WTNW stateR, FALSE
#define YES_WTNW_EMPTY stateR, TRUE
#define YES_WTMODEL mR, FALSE
#define YES_WTMODEL_NOINIT_S mR, TRUE
#define YES_WTMHPROPOSAL pR
#define YES_WTLASTTOGGLE length(time)>0, asInteger(time), INTEGER(lasttoggle)

#define NO_WTMODEL NULL, FALSE
#define NO_WTMHPROPOSAL  NULL
#define NO_WTLASTTOGGLE FALSE, 0, NULL

typedef struct{
  double *stats;
  WtNetwork *nwp;
  WtModel *m;
  WtMHProposal *MHp;
} WtErgmState;

WtErgmState *WtErgmStateInit(// Network settings
			     SEXP stateR, Rboolean empty,
                             // Model settings
			     SEXP mR, Rboolean noinit_s,
                             // Proposal settings
                             SEXP pR,
                             // Network state
                             Rboolean timings, int time, int *lasttoggle);
SEXP WtErgmStateRSave(SEXP startR, WtErgmState *s);
void WtErgmStateDestroy(WtErgmState *s);

#endif
