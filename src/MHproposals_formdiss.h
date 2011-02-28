#ifndef MHproposals_formdiss_H
#define MHproposals_formdiss_H

#include "edgetree.h"
#include "changestat.h"
#include "model.h"
#include "MHproposal.h"

void MH_Formation(MHproposal *MHp, Network *nwp);
void MH_FormationMLE(MHproposal *MHp, Network *nwp);
void MH_FormationTNT(MHproposal *MHp, Network *nwp);
void MH_FormationNonObservedMLE(MHproposal *MHp, Network *nwp);
/*void MH_DissolutionTNT(MHproposal *MHp, Network *nwp); */
void MH_Dissolution(MHproposal *MHp, Network *nwp);
void MH_DissolutionMLE(MHproposal *MHp, Network *nwp);
void MH_DissolutionNonObservedMLE(MHproposal *MHp, Network *nwp);

void MH_BipartiteFormation (MHproposal *MHp, Network *nwp);
void MH_BipartiteFormationTNT (MHproposal *MHp, Network *nwp);
void MH_BipartiteDissolution (MHproposal *MHp, Network *nwp);

#endif 

