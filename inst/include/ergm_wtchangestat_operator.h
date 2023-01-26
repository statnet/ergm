/*  File inst/include/ergm_wtchangestat_operator.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2023 Statnet Commons
 */
#ifndef _ERGM_WTCHANGESTAT_OPERATOR_H_
#define _ERGM_WTCHANGESTAT_OPERATOR_H_

#include "ergm_wtmodel.h"

#define WtSELECT_C_OR_D_BASED_ON_SUBMODEL(m)            \
  {                                                     \
    Rboolean c = TRUE;                                  \
    WtFOR_EACH_TERM(m) if(mtp->d_func) c=FALSE;         \
    if(c) mtp->d_func=NULL; else mtp->c_func=NULL;      \
  }

#define WtDELETE_IF_UNUSED_IN_SUBMODEL(what, m)                 \
  {                                                             \
    Rboolean used = FALSE;                                      \
    WtFOR_EACH_TERM((WtModel*)m) if(mtp->what) used=TRUE;       \
    if(!used) mtp->what=NULL;                                   \
  }

#define WtDELETE_IF_UNUSED_IN_SUBMODELS(what, ms, n)                    \
  {                                                                     \
    Rboolean used = FALSE;                                              \
    for(unsigned int i=0; i<n; i++)                                     \
      if(((WtModel**)ms)[i]) WtFOR_EACH_TERM(((WtModel**)ms)[i]) if(mtp->what) used=TRUE; \
    if(!used) mtp->what=NULL;                                           \
  }

#define WtPROPAGATE_X_SIGNAL(onwp, m) WtSEND_X_SIGNAL_INTO((onwp), (m), NULL, (m)->workspace, type, data);

#define WtPROPAGATE_X_SIGNAL_INTO(onwp, m, output)                      \
  {                                                                     \
    memset(output, 0, (m)->n_stats*sizeof(double));                     \
    WtSEND_X_SIGNAL_INTO((onwp), (m), NULL, (output), type, data);      \
  }

#define WtPROPAGATE_X_SIGNAL_ADDONTO(onwp, m, output)                   \
  {                                                                     \
    WtPROPAGATE_X_SIGNAL_INTO((onwp), (m), (m)->workspace);             \
    addonto((output), (m)->workspace, (m)->n_stats);                    \
  }

#define WtX_CHANGESTAT_PROPAGATE_FN(a, getstorage, getm)                \
  WtX_CHANGESTAT_FN(a) {                                                \
    getstorage;                                                         \
    WtModelTerm *_mymtp = mtp;                                          \
    WtSEND_X_SIGNAL_INTO(nwp, (getm), NULL, _mymtp->dstats, type, data); \
  }

#endif // _ERGM_WTCHANGESTAT_OPERATOR_H_
