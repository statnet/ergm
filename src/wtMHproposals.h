/*  File src/wtMHproposals.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2018 Statnet Commons
 */
#ifndef WTMHPROPOSALS_H
#define WTMHPROPOSALS_H

#include "wtMHproposal.h"

void MH_Unif(WtMHproposal *MHp, WtNetwork *nwp);
void MH_UnifNonObserved(WtMHproposal *MHp, WtNetwork *nwp);
void MH_DiscUnif(WtMHproposal *MHp, WtNetwork *nwp);
void MH_DiscUnifNonObserved(WtMHproposal *MHp, WtNetwork *nwp);

#endif 



