/*  File inst/include/inc/ergm_state.h.template.do_not_include_directly.h in
 *  package ergm, part of the Statnet suite of packages for network analysis,
 *  https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */

#include "../ergm_constants.h"

typedef struct{
  SEXP R;
  double *stats;
  ETYPE(Network) *nwp;
  ETYPE(Model) *m;
  ETYPE(MHProposal) *MHp;
  SEXP save;
} ETYPE(ErgmState);

ETYPE(ErgmState) *ETYPE(ErgmStateInit)(SEXP stateR,
                             unsigned int flags);
SEXP ETYPE(ErgmStateRSave)(ETYPE(ErgmState) *s);
void ETYPE(ErgmStateDestroy)(ETYPE(ErgmState) *s);
SEXP ETYPE(ErgmStateArrayClear)(void);
