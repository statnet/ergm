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

#endif
