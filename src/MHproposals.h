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


/* Maximum tries (up to an MH-specific constant). */
#define MAX_TRIES 5000


/* MH_* proposal failed codes. */
/* Heads: */
#define MH_FAILED 0
/* Tails: */
#define MH_UNRECOVERABLE 0
#define MH_IMPOSSIBLE 1
#define MH_UNSUCCESSFUL 2


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

void MH_Formation(MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_FormationTNT(MHproposal *MHp, DegreeBound *bd, Network *nwp);
//void MH_DissolutionTNT(MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_Dissolution(MHproposal *MHp, DegreeBound *bd, Network *nwp);
#endif 


/* Summary of which functions can currently be called from ergm (Oct. 9 2007)
MH_randomtoggle : Done : 
   ergm (..., constraint = "none", control=ergm.control(prop.weights="random"))

MH_TNT : Done :
   ergm (..., constraint = "none", control=ergm.control(prop.weights="default"))
or ergm (..., constraint = "none", control=ergm.control(prop.weights="TNT"))

MH_ConstantEdges : Done :
   ergm (..., constraint = "edges", control=ergm.control(prop.weights="default"))
or ergm (..., constraint = "edges", control=ergm.control(prop.weights="random"))

MH_CondDegTetra : Not Done
MH_CondDegreeDist : Not Done *** What is the difference with CondDegree? ***
MH_CondOutDegreeDist : Done :
   ergm (..., constraint = "outdegrees", control=ergm.control(prop.weights="default"))
or ergm (..., constraint = "outdegrees", control=ergm.control(prop.weights="random"))
MH_CondInDegreeDist : Done :
   ergm (..., constraint = "indegrees", control=ergm.control(prop.weights="default"))
or ergm (..., constraint = "indegrees", control=ergm.control(prop.weights="random"))

MH_CondDegree : Done : *** Should it work for directed networks? ***
   ergm (..., constraint = "degrees", control=ergm.control(prop.weights="default"))
or ergm (..., constraint = "degrees", control=ergm.control(prop.weights="random"))
MH_CondDegHexadToggles : Done : called from CondDegree with prob 0.1
MH_CondDegTetradToggles : Done : called from CondDegree with prob 0.9

*** Note:  We are currently missing the following functions in C: 
    (a) MH_Hamming 
    (b) MH_HammingConstantEdges function in C, 
    (c) MH_nobetweengroupties
    though one can currently invoke (a) and (b) using 
   ergm (..., constraint = "hamming", control=ergm.control(prop.weights="default"))
or ergm (..., constraint = "hamming", control=ergm.control(prop.weights="random"))
or ergm (..., constraint = "edges+hamming", control=ergm.control(prop.weights="default"))
or ergm (..., constraint = "edges+hamming", control=ergm.control(prop.weights="random"))
but for (c), only InitMHP.nobetweengroupties exists

MH_RandomNode 
MH_randomtoggleNonObserved : Not Done :


MH_ConstrainedCondOutDegDist 
MH_OneRandomTnTNode 
MH_TwoRandomToggles 
MH_NodePairedTiesToggles 
MH_AllTogglesForOneNode 
MH_ReallocateWithReplacement 
MH_SwitchLabelTwoNodesToggles 
MH_ConstantEdgesToggles 
MH_OneConstrainedRandomToggle 
MH_ConstrainedTwoRandomToggles 
MH_ConstrainedNodePairedTiesToggles 
MH_ConstrainedCondDegDist 
MH_ConstrainedCondDegSwitchToggles 
MH_ConstrainedAllTogglesForOneNode 
MH_ConstrainedReallocateWithReplacement 
MH_ConstrainedSwitchLabelTwoNodesToggles 

MH_Formation
   ergm (..., constraint = "none", control=ergm.control(prop.weights="random"))
MH_FormationTNT : Done :                                                                              
   ergm (..., constraint = "none", control=ergm.control(prop.weights="default"))
or ergm (..., constraint = "none", control=ergm.control(prop.weights="TNT"))
//MH_DissolutionTNT
MH_Dissolution : Done :
   ergm (..., constraint = "none", control=ergm.control(prop.weights="default"))
or ergm (..., constraint = "none", control=ergm.control(prop.weights="TNT"))

------------------------------------------------
From MHproposals_bipartite.h:

MH_Bipartiterandomtoggle           
MH_BipartiteConstantEdges

MH_BipartiteHammingConstantEdges : Done (assuming network is bipartite):
   ergm (..., constraint = "edges+hamming", control=ergm.control(prop.weights="default"))
or ergm (..., constraint = "edges+hamming", control=ergm.control(prop.weights="random"))

MH_BipartiteHamming : Done (assuming network is bipartite):
   ergm (..., constraint = "edges+hamming", control=ergm.control(prop.weights="default"))
or ergm (..., constraint = "edges+hamming", control=ergm.control(prop.weights="random"))

MH_BipartiteCondDegreeDist         

MH_BipartiteFormation : Done (assuming network is bipartite):
   ergm (..., constraint = "none", control=ergm.control(prop.weights="random"))
MH_BipartiteFormationTNT : Done (assuming network is bipartite): 
   ergm (..., constraint = "none", control=ergm.control(prop.weights="default"))
or ergm (..., constraint = "none", control=ergm.control(prop.weights="TNT"))

MH_BipartiterandomtoggleNonObserved : Not Done (assuming network is bipartite):

MH_BipartiteDissolution : Done (assuming network is bipartite):
   ergm (..., constraint = "none", control=ergm.control(prop.weights="default"))
or ergm (..., constraint = "none", control=ergm.control(prop.weights="TNT"))

*/

