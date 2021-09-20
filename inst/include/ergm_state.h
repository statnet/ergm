/*  File inst/include/ergm_state.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2021 Statnet Commons
 */
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
  SEXP save;
} ErgmState;

ErgmState *ErgmStateInit(// Network settings
                         SEXP stateR,
                         unsigned int flags);
SEXP ErgmStateRSave(ErgmState *s);
void ErgmStateDestroy(ErgmState *s);
SEXP ErgmStateArrayClear();

#endif
