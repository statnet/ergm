/*  File inst/include/ergm_changestat_operator.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */
#ifndef _ERGM_CHANGESTAT_OPERATOR_H_
#define _ERGM_CHANGESTAT_OPERATOR_H_

#include "ergm_model.h"

#define SELECT_C_OR_D_BASED_ON_SUBMODEL(m)              \
  {                                                     \
    Rboolean c = TRUE;                                  \
    FOR_EACH_TERM((Model*)m) if(mtp->d_func) c=FALSE;   \
    if(c) mtp->d_func=NULL; else mtp->c_func=NULL;      \
  }

#define DELETE_IF_UNUSED_IN_SUBMODEL(what, m)           \
  {                                                     \
    Rboolean used = FALSE;                              \
    FOR_EACH_TERM((Model*)m) if(mtp->what) used=TRUE;   \
    if(!used) mtp->what=NULL;                           \
  }

#define DELETE_IF_UNUSED_IN_SUBMODELS(what, ms, n)                      \
  {                                                                     \
    Rboolean used = FALSE;                                              \
    for(unsigned int i=0; i<n; i++)                                     \
      if(((Model**)ms)[i]) FOR_EACH_TERM(((Model**)ms)[i]) if(mtp->what) used=TRUE; \
    if(!used) mtp->what=NULL;                                           \
  }

#define PROPAGATE_X_SIGNAL(onwp, m) SEND_X_SIGNAL_INTO((onwp), (m), NULL, (m)->workspace, type, data);

#define PROPAGATE_X_SIGNAL_INTO(onwp, m, output)                        \
  {                                                                     \
    memset(output, 0, (m)->n_stats*sizeof(double));                     \
    SEND_X_SIGNAL_INTO((onwp), (m), NULL, (output), type, data);        \
  }

#define PROPAGATE_X_SIGNAL_ADDONTO(onwp, m, output)                     \
  {                                                                     \
    PROPAGATE_X_SIGNAL_INTO((onwp), (m), (m)->workspace);               \
    addonto((output), (m)->workspace, (m)->n_stats);                    \
  }

#define X_CHANGESTAT_PROPAGATE_FN(a, getstorage, getm)                  \
  X_CHANGESTAT_FN(a) {                                                  \
    getstorage;                                                         \
    ModelTerm *_mymtp = mtp;                                            \
    SEND_X_SIGNAL_INTO(nwp, (getm), NULL, _mymtp->dstats, type, data);   \
  }

#endif // _ERGM_CHANGESTAT_OPERATOR_H_
