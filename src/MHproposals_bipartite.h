#ifndef MHproposals_bipartite_H
#define MHproposals_bipartite_H

#include "edgetree.h"
#include "changestat.h"
#include "model.h"
#include "MHproposal.h"

void MH_BipartiteHammingConstantEdges (MHproposal *MHp, Network *nwp);
void MH_BipartiteHammingTNT (MHproposal *MHp, Network *nwp);
void MH_BipartiteCondDegreeDist (MHproposal *MHp, Network *nwp);

#endif 


