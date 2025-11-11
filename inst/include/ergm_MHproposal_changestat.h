#ifndef _ERGM_MHPROPOSAL_CHANGESTAT_H_
#define _ERGM_MHPROPOSAL_CHANGESTAT_H_

#include "ergm_model.h"

#include "ergm_edgetype_set_binary.h"

#include "inc/ergm_MHproposal_changestat.h.template.do_not_include_directly.h"

#include "ergm_edgetype_unset.h"

#define _CHECK_CHANGESTATS2(m, edgestate)                       \
  if (m && ChangeStats1_changed(MHp, nwp, m, edgestate)) {      \
    MHp->toggletail[0]=MH_FAILED;                               \
    MHp->togglehead[0]=MH_CONSTRAINT;                           \
    return;                                                     \
  }

#define _CHECK_CHANGESTATS1(m)                  \
  if (m && ChangeStats_changed(MHp, nwp, m)) {  \
    MHp->toggletail[0]=MH_FAILED;               \
    MHp->togglehead[0]=MH_CONSTRAINT;           \
    return;                                     \
  }

#define CHECK_CHANGESTATS(...) _GET_OVERRIDE2(__VA_ARGS__, _CHECK_CHANGESTATS2, _CHECK_CHANGESTATS1,)(__VA_ARGS__)

#endif
