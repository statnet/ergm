#ifndef _ERGM_STATE_H_
#define _ERGM_STATE_H_

#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_MHproposal.h"
#include "ergm_model.h"

#define ARGS_STATE SEXP stateR
#define ARGS_LASTTOGGLE SEXP time, SEXP lasttoggle

#define YES_STATE stateR, FALSE, FALSE
#define YES_STATE_EMPTY_NO_INIT_S stateR, TRUE, TRUE
#define YES_LASTTOGGLE length(time)>0, asInteger(time), INTEGER(lasttoggle)

#define NO_LASTTOGGLE FALSE, 0, NULL

typedef struct{
  double *stats;
  Network *nwp;
  Model *m;
  MHProposal *MHp;
} ErgmState;

ErgmState *ErgmStateInit(// Network settings
                         SEXP stateR,
                         Rboolean empty, Rboolean noinit_s,
                         // Network state
                         Rboolean timings, int time, int *lasttoggle);
SEXP ErgmStateRSave(SEXP startR, ErgmState *s);
void ErgmStateDestroy(ErgmState *s);

#endif
