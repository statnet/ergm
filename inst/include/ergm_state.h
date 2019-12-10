#ifndef _ERGM_STATE_H_
#define _ERGM_STATE_H_

#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_MHproposal.h"
#include "ergm_model.h"

#define ARGS_NW SEXP stateR
#define ARGS_MODEL SEXP mR
#define ARGS_MHPROPOSAL SEXP pR
#define ARGS_LASTTOGGLE SEXP time, SEXP lasttoggle

#define YES_NW stateR, FALSE
#define YES_NW_EMPTY stateR, TRUE
#define YES_MODEL mR, FALSE
#define YES_MODEL_NOINIT_S mR, TRUE
#define YES_MHPROPOSAL pR
#define YES_LASTTOGGLE length(time)>0, asInteger(time), INTEGER(lasttoggle)

#define NO_MODEL NULL, FALSE
#define NO_MHPROPOSAL  NULL
#define NO_LASTTOGGLE FALSE, 0, NULL

typedef struct{
  double *stats;
  Network *nwp;
  Model *m;
  MHProposal *MHp;
} ErgmState;

ErgmState *ErgmStateInit(// Network settings
                         SEXP stateR, Rboolean empty,
                         // Model settings
                         SEXP mR, Rboolean noinit_s,
                         // Proposal settings
                         SEXP pR,
                         // Network state
                         Rboolean timings, int time, int *lasttoggle);
SEXP ErgmStateRSave(SEXP startR, ErgmState *s);
void ErgmStateDestroy(ErgmState *s);

#endif
