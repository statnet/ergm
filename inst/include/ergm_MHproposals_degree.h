/*  File inst/include/ergm_MHproposals_degree.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */
#include "ergm_MHproposal.h"
#include "ergm_changestat.h"

#ifdef __cplusplus
extern "C" {
#endif

MH_P_FN(MH_CondDegreeTetrad);
MH_P_FN(MH_CondDegreeHexad);
MH_P_FN(MH_CondDegree);
MH_P_FN(MH_CondOutDegree);
MH_P_FN(MH_CondInDegree);
MH_P_FN(MH_CondB1Degree);
MH_P_FN(MH_CondB2Degree);

#ifdef __cplusplus
}
#endif
