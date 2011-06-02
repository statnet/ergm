#ifndef MHproposals_DynMLE_H
#define MHproposals_DynMLE_H

#include "edgetree.h"
#include "changestat.h"
#include "model.h"
#include "MHproposal.h"

void MH_FormationMLE(MHproposal *MHp, Network *nwp);
void MH_FormationNonObservedMLE(MHproposal *MHp, Network *nwp);
void MH_DissolutionMLE(MHproposal *MHp, Network *nwp);
void MH_DissolutionNonObservedMLE(MHproposal *MHp, Network *nwp);

#endif 

