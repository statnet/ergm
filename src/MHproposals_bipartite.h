/*  File src/MHproposals_bipartite.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2018 Statnet Commons
 */
#ifndef MHproposals_bipartite_H
#define MHproposals_bipartite_H

#include "edgetree.h"
#include "changestat.h"
#include "model.h"
#include "MHproposal.h"

void MH_BipartiteHammingConstantEdges (MHproposal *MHp, Network *nwp);
void MH_BipartiteHammingTNT (MHproposal *MHp, Network *nwp);
void MH_BipartiteCondDegreeDist (MHproposal *MHp, Network *nwp);

#endif 


