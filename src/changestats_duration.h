#ifndef CHANGESTATS_duration_H
#define CHANGESTATS_duration_H

#include "edgeTree.h"
#include "basechangeStats.h"

void d_D_on (int ntoggles, Vertex *heads, Vertex *tail, 
                ModelTerm *mtp, Network *nwp);
void d_D_off (int ntoggles, Vertex *heads, Vertex *tail, 
                ModelTerm *mtp, Network *nwp);
double mean_duration(Network *nwp);

#endif








