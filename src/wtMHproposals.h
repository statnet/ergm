#ifndef WTMHPROPOSALS_H
#define WTMHPROPOSALS_H

#include "wtMHproposal.h"

void MH_Unif(WtMHproposal *MHp, WtNetwork *nwp);
void MH_UnifNonObserved(WtMHproposal *MHp, WtNetwork *nwp);
void MH_DiscUnif(WtMHproposal *MHp, WtNetwork *nwp);
void MH_DiscUnifNonObserved(WtMHproposal *MHp, WtNetwork *nwp);

#endif 



