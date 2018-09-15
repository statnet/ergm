/*  File src/MHProposals.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2017 Statnet Commons
 */
#ifndef MHProposals_H
#define MHProposals_H

#include "MHproposal.h"

void MH_randomtoggle (MHProposal *MHp, Network *nwp);
void MH_TNT (MHProposal *MHp, Network *nwp);
void MH_TNT10 (MHProposal *MHp, Network *nwp);
void MH_ConstantEdges (MHProposal *MHp, Network *nwp);
void MH_CondDegreeDist (MHProposal *MHp, Network *nwp);
void MH_CondOutDegreeDist (MHProposal *MHp, Network *nwp);
void MH_CondInDegreeDist (MHProposal *MHp, Network *nwp);
void MH_RandomNode (MHProposal *MHp, Network *nwp);
void MH_randomtoggleList (MHProposal *MHp, Network *nwp);

void MH_ConstrainedCondOutDegDist (MHProposal *MHp, Network *nwp);
void MH_OneRandomTnTNode (MHProposal *MHp, Network *nwp);
void MH_TwoRandomToggles (MHProposal *MHp, Network *nwp);
void MH_NodePairedTiesToggles (MHProposal *MHp, Network *nwp);
void MH_AllTogglesForOneNode (MHProposal *MHp, Network *nwp);
void MH_ReallocateWithReplacement (MHProposal *MHp, Network *nwp);
void MH_SwitchLabelTwoNodesToggles (MHProposal *MHp, Network *nwp);
void MH_ConstantEdgesToggles (MHProposal *MHp, Network *nwp);
void MH_OneConstrainedRandomToggle (MHProposal *MHp, Network *nwp);
void MH_ConstrainedTwoRandomToggles (MHProposal *MHp, Network *nwp);
void MH_ConstrainedNodePairedTiesToggles (MHProposal *MHp, Network *nwp);
void MH_ConstrainedCondDegDist (MHProposal *MHp, Network *nwp);
void MH_ConstrainedCondDegSwitchToggles (MHProposal *MHp, Network *nwp);
void MH_ConstrainedAllTogglesForOneNode (MHProposal *MHp, Network *nwp);
void MH_ConstrainedReallocateWithReplacement (MHProposal *MHp, Network *nwp);
void MH_ConstrainedSwitchLabelTwoNodesToggles (MHProposal *MHp, Network *nwp);

#endif 



