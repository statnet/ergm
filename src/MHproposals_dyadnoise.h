/*  File src/MHproposals_dyadnoise.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2017 Statnet Commons
 */
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


