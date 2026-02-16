/*  File inst/include/ergm_wtMHproposal_changestat.h in package ergm, part of
 *  the Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */
#ifndef _ERGM_WTMHPROPOSAL_CHANGESTAT_H_
#define _ERGM_WTMHPROPOSAL_CHANGESTAT_H_

#include "ergm_wtmodel.h"

#include "ergm_edgetype_set_double.h"

#include "inc/ergm_MHproposal_changestat.h.template.do_not_include_directly.h"

#include "ergm_edgetype_unset.h"


#define _WtCHECK_CHANGESTATS2(m, edgestate) {                   \
  if (m && WtChangeStats1_changed(MHp, nwp, m, edgestate)) {    \
    MHp->toggletail[0]=MH_FAILED;                               \
    MHp->togglehead[0]=MH_CONSTRAINT;                           \
    return;                                                     \
  }

#define _WtCHECK_CHANGESTATS1(m) {                      \
  if (m && WtChangeStats_changed(MHp, nwp, m)) {        \
    MHp->toggletail[0]=MH_FAILED;                       \
    MHp->togglehead[0]=MH_CONSTRAINT;                   \
    return;                                             \
  }

#define WtCHECK_CHANGESTATS(...) _GET_OVERRIDE2(__VA_ARGS__, _WtCHECK_CHANGESTATS2, _WtCHECK_CHANGESTATS1,)(__VA_ARGS__)

#ifdef STUBFILE
#define STRICT_Wt_HEADERS
#endif

#ifndef STRICT_Wt_HEADERS
#define CHECK_CHANGESTATS WtCHECK_CHANGESTATS
#endif // STRICT_Wt_HEADERS

#endif // _ERGM_WTMHPROPOSAL_CHANGESTAT_H_
