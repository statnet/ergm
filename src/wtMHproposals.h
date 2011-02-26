#ifndef WTMHPROPOSALS_H
#define WTMHPROPOSALS_H

#include "wtMHproposal.h"

void MH_Poisson(WtMHproposal *MHp, WtNetwork *nwp);
void MH_PoissonNonObserved(WtMHproposal *MHp, WtNetwork *nwp);
void MH_CompleteOrdering(WtMHproposal *MHp, WtNetwork *nwp);

#endif 



