#ifndef WTMHPROPOSALS_BIPARTITE_H
#define WTMHPROPOSALS_BIPARTITE_H

#include "wtMHproposal.h"

void MH_BipartitePoisson(WtMHproposal *MHp, WtNetwork *nwp);
void MH_BipartitePoissonNonObserved(WtMHproposal *MHp, WtNetwork *nwp);
void MH_CompleteOrderingBipartite(WtMHproposal *MHp, WtNetwork *nwp);

#endif 

