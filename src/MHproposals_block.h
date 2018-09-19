/*  File src/MHProposals_block.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2017 Statnet Commons
 */
#ifndef MHProposals_Block_H
#define MHProposals_Block_H

#include "ergm_MHproposal.h"

void MH_blockdiag (MHProposal *MHp, Network *nwp);
void MH_blockdiagTNT (MHProposal *MHp, Network *nwp);

#endif 



