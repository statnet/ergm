/*  File src/MHproposals.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2020 Statnet Commons
 */
#ifndef MHProposals_H
#define MHProposals_H

#include "ergm_MHproposal.h"

MH_P_FN(MH_randomtoggle);
MH_P_FN(MH_TNT);
MH_P_FN(MH_TNT10);
MH_P_FN(MH_ConstantEdges);
MH_P_FN(MH_CondDegreeDist);
MH_P_FN(MH_CondOutDegreeDist);
MH_P_FN(MH_CondInDegreeDist);
MH_P_FN(MH_RandomNode);
MH_P_FN(MH_randomtoggleList);

MH_P_FN(MH_ConstrainedCondOutDegDist);
MH_P_FN(MH_OneRandomTnTNode);
MH_P_FN(MH_TwoRandomToggles);
MH_P_FN(MH_NodePairedTiesToggles);
MH_P_FN(MH_AllTogglesForOneNode);
MH_P_FN(MH_ReallocateWithReplacement);
MH_P_FN(MH_SwitchLabelTwoNodesToggles);
MH_P_FN(MH_ConstantEdgesToggles);
MH_P_FN(MH_OneConstrainedRandomToggle);
MH_P_FN(MH_ConstrainedTwoRandomToggles);
MH_P_FN(MH_ConstrainedNodePairedTiesToggles);
MH_P_FN(MH_ConstrainedCondDegDist);
MH_P_FN(MH_ConstrainedCondDegSwitchToggles);
MH_P_FN(MH_ConstrainedAllTogglesForOneNode);
MH_P_FN(MH_ConstrainedReallocateWithReplacement);
MH_P_FN(MH_ConstrainedSwitchLabelTwoNodesToggles);

#endif 



