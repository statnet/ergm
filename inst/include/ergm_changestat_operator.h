#ifndef _ERGM_CHANGESTAT_OPERATOR_H_
#define _ERGM_CHANGESTAT_OPERATOR_H_

#include "ergm_model.h"

#define SELECT_C_OR_D_BASED_ON_SUBMODEL(m)              \
  {                                                     \
    Rboolean c = TRUE;                                  \
    FOR_EACH_TERM(m) if(mtp->d_func) c=FALSE;           \
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

#endif // _ERGM_CHANGESTAT_OPERATOR_H_
