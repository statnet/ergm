/*  File src/MHProposals_dyadnoise.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2017 Statnet Commons
 */
#ifndef MHProposals_dyadnoise_H
#define MHProposals_dyadnoise_H

#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_model.h"
#include "ergm_MHproposal.h"

void MH_dyadnoiseTNT (MHProposal *MHp, Network *nwp);
void MH_dyadnoisemTNT (MHProposal *MHp, Network *nwp);
void MH_dyadnoise (MHProposal *MHp, Network *nwp);
void MH_dyadnoisem (MHProposal *MHp, Network *nwp);

#endif 


