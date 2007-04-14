#ifndef MHproposals_H
#define MHproposals_H

#include "edgeTree.h"
#include "basechangeStats.h"
#include "model.h"
#include "MCMC.h"

#define NO_EDGE       0x00 /*these four used in realocateWithReplacement */
#define OLD_EDGE      0x01 
#define NEW_EDGE      0x02 
#define CAN_IGNORE    (OLD_EDGE | NEW_EDGE)  

void MH_randomtoggle (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_TNT (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_ConstantEdges (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_CondDegTetra (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_CondDegreeDist (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_CondOutDegreeDist (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_CondInDegreeDist (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_CondDegree (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_CondDegHexadToggles (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_CondDegTetradToggles (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_RandomNode (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_randomtoggleNotObserved (MHproposal *MHp, DegreeBound *bd, Network *nwp);

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

