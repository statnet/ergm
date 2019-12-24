#ifndef _ERGM_CHANGESTAT_OPERATOR_H_
#define _ERGM_CHANGESTAT_OPERATOR_H_

#include "ergm_model.h"

#define SELECT_C_OR_D_BASED_ON_SUBMODEL(m)              \
  {                                                     \
    Rboolean c = TRUE;                                  \
    FOR_EACH_TERM(m) if(mtp->d_func) c=FALSE;           \
    if(c) mtp->d_func=NULL; else mtp->c_func=NULL;      \
  }

#endif // _ERGM_CHANGESTAT_OPERATOR_H_
