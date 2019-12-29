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

#endif // _ERGM_WTCHANGESTAT_OPERATOR_H_
