#ifndef _ERGM_STATE_H_
#define _ERGM_STATE_H_

#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_MHproposal.h"
#include "ergm_model.h"

#define ARGS_STATE SEXP stateR

#define YES_STATE stateR, FALSE, FALSE
#define YES_STATE_EMPTY_NO_INIT_S stateR, TRUE, TRUE

typedef struct{
  double *stats;
  Network *nwp;
  Model *m;
  MHProposal *MHp;
} ErgmState;

ErgmState *ErgmStateInit(// Network settings
                         SEXP stateR,
                         Rboolean empty, Rboolean noinit_s);
SEXP ErgmStateRSave(SEXP startR, ErgmState *s);
void ErgmStateDestroy(ErgmState *s);

#endif
