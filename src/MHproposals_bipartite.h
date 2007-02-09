#ifndef MHproposals_bipartite_H
#define MHproposals_bipartite_H

#include "edgeTree.h"
#include "basechangeStats.h"
#include "model.h"
#include "MCMC.h"

void MH_Bipartiterandomtoggle (MHproposal *MHp,  DegreeBound *bd, Network *nwp);
void MH_BipartiteConstantEdges (MHproposal *MHp,  DegreeBound *bd, Network *nwp);
void MH_BipartiteHammingConstantEdges (MHproposal *MHp,  DegreeBound *bd, Network *nwp);
void MH_BipartiteHamming (MHproposal *MHp,  DegreeBound *bd, Network *nwp);
void MH_BipartiteCondDegreeDist (MHproposal *MHp, DegreeBound *bd, Network *nwp);
void MH_BipartiteFormation (MHproposal *MHp,  DegreeBound *bd, Network *nwp);
void MH_BipartiteFormationTNT (MHproposal *MHp,  DegreeBound *bd, Network *nwp);
 
#endif 


