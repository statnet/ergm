#ifndef _ERGM_STATE_H_
#define _ERGM_STATE_H_

#include "ergm_constants.h"
#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_MHproposal.h"
#include "ergm_model.h"

typedef struct{
  SEXP R;
  double *stats;
  Network *nwp;
  Model *m;
  MHProposal *MHp;
} ErgmState;

ErgmState *ErgmStateInit(// Network settings
                         SEXP stateR,
                         unsigned int flags);
SEXP ErgmStateRSave(ErgmState *s);
void ErgmStateDestroy(ErgmState *s);
SEXP ErgmStateArrayClear();

#endif
