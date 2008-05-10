#ifndef CHANGESTATS_DURATION_H
#define CHANGESTATS_DURATION_H

#include "wtedgetree.h"
#include "changestats.h"

void d_D_on (int ntoggles, Vertex *heads, Vertex *tail, 
                ModelTerm *mtp, Network *nwp);
void d_D_off (int ntoggles, Vertex *heads, Vertex *tail, 
                ModelTerm *mtp, Network *nwp);
double mean_duration(Network *nwp);
D_CHANGESTAT_FN(d_edges_ageinterval);
#endif








