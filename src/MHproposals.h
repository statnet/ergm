/*  File src/MHproposals.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2017 Statnet Commons
 */
#ifndef MHproposals_H
#define MHproposals_H

#include "MHproposal.h"

void MH_randomtoggle (MHproposal *MHp, Network *nwp);
void MH_TNT (MHproposal *MHp, Network *nwp);
void MH_TNT10 (MHproposal *MHp, Network *nwp);
void MH_ConstantEdges (MHproposal *MHp, Network *nwp);
void MH_CondDegreeDist (MHproposal *MHp, Network *nwp);
void MH_CondOutDegreeDist (MHproposal *MHp, Network *nwp);
void MH_CondInDegreeDist (MHproposal *MHp, Network *nwp);
void MH_RandomNode (MHproposal *MHp, Network *nwp);
void MH_randomtoggleList (MHproposal *MHp, Network *nwp);

void MH_ConstrainedCondOutDegDist (MHproposal *MHp, Network *nwp);
void MH_OneRandomTnTNode (MHproposal *MHp, Network *nwp);
void MH_TwoRandomToggles (MHproposal *MHp, Network *nwp);
void MH_NodePairedTiesToggles (MHproposal *MHp, Network *nwp);
void MH_AllTogglesForOneNode (MHproposal *MHp, Network *nwp);
void MH_ReallocateWithReplacement (MHproposal *MHp, Network *nwp);
void MH_SwitchLabelTwoNodesToggles (MHproposal *MHp, Network *nwp);
void MH_ConstantEdgesToggles (MHproposal *MHp, Network *nwp);
void MH_OneConstrainedRandomToggle (MHproposal *MHp, Network *nwp);
void MH_ConstrainedTwoRandomToggles (MHproposal *MHp, Network *nwp);
void MH_ConstrainedNodePairedTiesToggles (MHproposal *MHp, Network *nwp);
void MH_ConstrainedCondDegDist (MHproposal *MHp, Network *nwp);
void MH_ConstrainedCondDegSwitchToggles (MHproposal *MHp, Network *nwp);
void MH_ConstrainedAllTogglesForOneNode (MHproposal *MHp, Network *nwp);
void MH_ConstrainedReallocateWithReplacement (MHproposal *MHp, Network *nwp);
void MH_ConstrainedSwitchLabelTwoNodesToggles (MHproposal *MHp, Network *nwp);

#endif 



