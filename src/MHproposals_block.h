/*  File src/MHproposals_block.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2018 Statnet Commons
 */
#ifndef MHproposals_Block_H
#define MHproposals_Block_H

#include "MHproposal.h"

void MH_blockdiag (MHproposal *MHp, Network *nwp);
void MH_blockdiagTNT (MHproposal *MHp, Network *nwp);

#endif 



