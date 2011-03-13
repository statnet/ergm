#ifndef WTCHANGESTATS_RANK_H
#define WTCHANGESTATS_RANK_H

#include "wtedgetree.h"
#include "wtchangestat.h"

#define IF_1_EGO_SWAPS_2_ALTERS(subroutine){if(ntoggles==2 && tails[0]==tails[1]){ZERO_ALL_CHANGESTATS(); Vertex t=tails[0], h1=heads[0], h2=heads[1]; {subroutine};}else{d_from_s(ntoggles, tails, heads, weights, mtp, nwp);};}

WtD_CHANGESTAT_FN(d_inconsistency_rank); WtS_CHANGESTAT_FN(s_inconsistency_rank);

WtD_CHANGESTAT_FN(d_deference); WtS_CHANGESTAT_FN(s_deference);

WtD_CHANGESTAT_FN(d_nodeicov_rank); WtS_CHANGESTAT_FN(s_nodeicov_rank);

WtD_CHANGESTAT_FN(d_nonconformity); WtS_CHANGESTAT_FN(s_nonconformity);

WtD_CHANGESTAT_FN(d_nonconformity_decay); WtS_CHANGESTAT_FN(s_nonconformity_thresholds);

WtD_CHANGESTAT_FN(d_nonconformity_thresholds); WtS_CHANGESTAT_FN(s_nonconformity_thresholds);

#endif
