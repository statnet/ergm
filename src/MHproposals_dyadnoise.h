#ifndef MHproposals_dyadnoise_H
#define MHproposals_dyadnoise_H

#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_model.h"
#include "ergm_MHproposal.h"

void MH_dyadnoiseTNT (MHproposal *MHp, Network *nwp);
void MH_dyadnoisemTNT (MHproposal *MHp, Network *nwp);
void MH_dyadnoise (MHproposal *MHp, Network *nwp);
void MH_dyadnoisem (MHproposal *MHp, Network *nwp);

#endif 


