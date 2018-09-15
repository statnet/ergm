/*  File src/MHProposals_degree.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2017 Statnet Commons
 */
#ifndef MHPROPOSALS_DEGREE_H
#define MHPROPOSALS_DEGREE_H
#include "MHproposal.h"

void MH_CondDegree (MHProposal *MHp, Network *nwp);
void MH_CondDegreeTetrad (MHProposal *MHp, Network *nwp);
void MH_CondDegreeHexad (MHProposal *MHp, Network *nwp);
void MH_CondOutDegree(MHProposal *MHp, Network *nwp);
void MH_CondInDegree(MHProposal *MHp, Network *nwp); 
void MH_CondB1Degree(MHProposal *MHp, Network *nwp); 
void MH_CondB2Degree(MHProposal *MHp, Network *nwp);  
void MH_CondDegreeMix(MHProposal *MHp, Network *nwp);
void MH_CondDegreeTetradMixMore(MHProposal *MHp, Network *nwp);
void MH_CondDegreeTetradMixLess(MHProposal *MHp, Network *nwp);
void MH_CondDegreeMixChangeOrig(MHProposal *MHp, Network *nwp);

#endif
