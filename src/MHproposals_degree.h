/*
 *  File ergm/src/MHproposals_degree.h
 *  Part of the statnet package, http://statnet.org
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) in
 *    http://statnet.org/attribution
 *
 *  Copyright 2012 the statnet development team
 */
#ifndef MHPROPOSALS_DEGREE_H
#define MHPROPOSALS_DEGREE_H
#include "MHproposal.h"

void MH_CondDegreeSimple (MHproposal *MHp, Network *nwp);
void MH_CondDegreeSimpleTetrad (MHproposal *MHp, Network *nwp);
void MH_CondDegreeSimpleHexad (MHproposal *MHp, Network *nwp);
#endif
