#ifndef WTCHANGESTATS_RANK_H
#define WTCHANGESTATS_RANK_H

#include "wtedgetree.h"
#include "wtchangestat.h"

#define IF_1_EGO_SWAPS_2_ALTERS(subroutine){if(ntoggles==2 && tails[0]==tails[1] && GETWT(tails[0], heads[0])==weights[1] && GETWT(tails[1], heads[1])==weights[0]){ZERO_ALL_CHANGESTATS(); Vertex t=tails[0], h1=heads[0], h2=heads[1]; {subroutine};}else{D_FROM_S};}

// Case 1: One ego changes the value of one alter.
// Case 2: One ego swaps the values of two alters.
// Case 3: The alters are adjacent?
#define OPTIMAL_RANK_D(case1sub,case2sub){if(ntoggles==1){ZERO_ALL_CHANGESTATS(); Vertex t=tails[0], h=heads[0]; {case1sub};} else if(ntoggles==2 && tails[0]==tails[1] && GETWT(tails[0], heads[0])==weights[1] && GETWT(tails[1], heads[1])==weights[0]){ZERO_ALL_CHANGESTATS(); Vertex t=tails[0], h1=heads[0], h2=heads[1]; {case2sub};} else{D_FROM_S};}

WtD_CHANGESTAT_FN(d_edgecov_rank); WtS_CHANGESTAT_FN(s_edgecov_rank);

WtD_CHANGESTAT_FN(d_inconsistency_rank); WtS_CHANGESTAT_FN(s_inconsistency_rank);

WtD_CHANGESTAT_FN(d_inconsistency_cov_rank); WtS_CHANGESTAT_FN(s_inconsistency_cov_rank);

WtD_CHANGESTAT_FN(d_deference); WtS_CHANGESTAT_FN(s_deference);

WtD_CHANGESTAT_FN(d_nodeicov_rank); WtS_CHANGESTAT_FN(s_nodeicov_rank);

WtD_CHANGESTAT_FN(d_nonconformity); WtS_CHANGESTAT_FN(s_nonconformity);

WtD_CHANGESTAT_FN(d_local_nonconformity); WtS_CHANGESTAT_FN(s_local_nonconformity);

WtD_CHANGESTAT_FN(d_nonconformity_decay); WtS_CHANGESTAT_FN(s_nonconformity_thresholds);

WtD_CHANGESTAT_FN(d_nonconformity_thresholds); WtS_CHANGESTAT_FN(s_nonconformity_thresholds);

WtD_CHANGESTAT_FN(d_tiedranks); WtS_CHANGESTAT_FN(s_tiedranks);

#endif
