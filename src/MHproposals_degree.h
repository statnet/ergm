/*  File src/MHproposals_degree.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2013 Statnet Commons
 */
#ifndef MHPROPOSALS_DEGREE_H
#define MHPROPOSALS_DEGREE_H
#include "MHproposal.h"

void MH_CondDegree (MHproposal *MHp, Network *nwp);
void MH_CondDegreeTetrad (MHproposal *MHp, Network *nwp);
void MH_CondDegreeHexad (MHproposal *MHp, Network *nwp);
void MH_CondOutDegree(MHproposal *MHp, Network *nwp);
void MH_CondInDegree(MHproposal *MHp, Network *nwp); 
void MH_CondB1Degree(MHproposal *MHp, Network *nwp); 
void MH_CondB2Degree(MHproposal *MHp, Network *nwp);  
void MH_CondDegreeMix(MHproposal *MHp, Network *nwp);
void MH_CondDegreeTetradMixMore(MHproposal *MHp, Network *nwp);
void MH_CondDegreeTetradMixLess(MHproposal *MHp, Network *nwp);
void MH_CondDegreeMixChangeOrig(MHproposal *MHp, Network *nwp);

#endif
