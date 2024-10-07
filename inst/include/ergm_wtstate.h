/*  File inst/include/ergm_wtstate.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2024 Statnet Commons
 */
#ifndef _ERGM_WTSTATE_H_
#define _ERGM_WTSTATE_H_

#include "ergm_constants.h"
#include "ergm_wtedgetree.h"
#include "ergm_wtchangestat.h"
#include "ergm_wtMHproposal.h"
#include "ergm_wtmodel.h"

typedef struct{
  SEXP R;
  double *stats;
  WtNetwork *nwp;
  WtModel *m;
  WtMHProposal *MHp;
  SEXP save;
} WtErgmState;

WtErgmState *WtErgmStateInit(SEXP stateR,
                             unsigned int flags);
SEXP WtErgmStateRSave(WtErgmState *s);
void WtErgmStateDestroy(WtErgmState *s);
SEXP WtErgmStateArrayClear(void);

#endif
