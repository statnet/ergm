/*
 *  File ergm/src/MHproposals_DynMoME.h
 *  Part of the statnet package, http://statnet.org
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) in
 *    http://statnet.org/attribution
 *
 *  Copyright 2012 the statnet development team
 */
#ifndef MHproposals_DynMoME_H
#define MHproposals_DynMoME_H

#include "edgetree.h"
#include "changestat.h"
#include "model.h"
#include "MHproposal.h"

void MH_Formation(MHproposal *MHp, Network *nwp);
void MH_FormationTNT(MHproposal *MHp, Network *nwp);
void MH_Dissolution(MHproposal *MHp, Network *nwp);

void MH_BipartiteFormation (MHproposal *MHp, Network *nwp);
void MH_BipartiteFormationTNT (MHproposal *MHp, Network *nwp);
void MH_BipartiteDissolution (MHproposal *MHp, Network *nwp);

#endif 

