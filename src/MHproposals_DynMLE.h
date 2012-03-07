/*
 *  File ergm/src/MHproposals_DynMLE.h
 *  Part of the statnet package, http://statnet.org
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) in
 *    http://statnet.org/attribution
 *
 *  Copyright 2012 the statnet development team
 */
#ifndef MHproposals_DynMLE_H
#define MHproposals_DynMLE_H

#include "edgetree.h"
#include "changestat.h"
#include "model.h"
#include "MHproposal.h"

void MH_FormationMLE(MHproposal *MHp, Network *nwp);
void MH_DissolutionMLE(MHproposal *MHp, Network *nwp);

void MH_FormationMLETNT(MHproposal *MHp, Network *nwp);

#endif 

