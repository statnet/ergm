#ifndef CHANGESTATS_duration_H
#define CHANGESTATS_duration_H

#include "edgeTree.h"
#include "basechangeStats.h"

void d_adegree (int ntoggles, Vertex *heads, Vertex *tails,
                ModelTerm *mtp, Network *nwp);
void d_edegree (int ntoggles, Vertex *heads, Vertex *tails,
                ModelTerm *mtp, Network *nwp);
void d_bimix (int ntoggles, Vertex *heads, Vertex *tails,
              ModelTerm *mtp, Network *nwp);
void d_formation (int ntoggles, Vertex *heads, Vertex *tails,
              ModelTerm *mtp, Network *nwp);
void d_dissolve (int ntoggles, Vertex *heads, Vertex *tails,
              ModelTerm *mtp, Network *nwp);

#endif
