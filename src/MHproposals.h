#ifndef MHproposals_H
#define MHproposals_H

#include "MHproposal.h"

void MH_randomtoggle (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_TNT (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_TNT10 (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_ConstantEdges (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_CondDegreeTetrad (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_CondDegreeDist (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_CondOutDegreeDist (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_CondInDegreeDist (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_CondDegree (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_CondDegreeHexadToggles (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_CondDegreeTetradToggles (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_RandomNode (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_randomtoggleNonObserved (MHproposal *MHp, DegreeBound *bd, Network *nwp);

void MH_ConstrainedCondOutDegDist (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_OneRandomTnTNode (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_TwoRandomToggles (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_NodePairedTiesToggles (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_AllTogglesForOneNode (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_ReallocateWithReplacement (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_SwitchLabelTwoNodesToggles (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_ConstantEdgesToggles (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_OneConstrainedRandomToggle (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_ConstrainedTwoRandomToggles (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_ConstrainedNodePairedTiesToggles (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_ConstrainedCondDegDist (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_ConstrainedCondDegSwitchToggles (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_ConstrainedAllTogglesForOneNode (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_ConstrainedReallocateWithReplacement (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_ConstrainedSwitchLabelTwoNodesToggles (MHproposal *MHp, DegreeBound *bd, Network *nwp);

#endif 



