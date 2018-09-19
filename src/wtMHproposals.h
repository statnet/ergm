/*  File src/wtProposals.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2017 Statnet Commons
 */
#ifndef WTMHPROPOSALS_H
#define WTMHPROPOSALS_H

#include "ergm_wtMHproposal.h"

void MH_Unif(WtMHProposal *MHp, WtNetwork *nwp);
void MH_UnifNonObserved(WtMHProposal *MHp, WtNetwork *nwp);
void MH_DiscUnif(WtMHProposal *MHp, WtNetwork *nwp);
void MH_DiscUnifNonObserved(WtMHProposal *MHp, WtNetwork *nwp);

#endif 



