#ifndef WTCHANGESTATS_RANK_H
#define WTCHANGESTATS_RANK_H

#include "wtedgetree.h"
#include "wtchangestat.h"

WtD_CHANGESTAT_FN(d_inconsistency_rank); WtS_CHANGESTAT_FN(s_inconsistency_rank);

WtD_CHANGESTAT_FN(d_deference); WtS_CHANGESTAT_FN(s_deference);

WtD_CHANGESTAT_FN(d_nodeicov_rank); WtS_CHANGESTAT_FN(s_nodeicov_rank);

WtD_CHANGESTAT_FN(d_nonconformity); WtS_CHANGESTAT_FN(s_nonconformity);

WtD_CHANGESTAT_FN(d_nonconformity_decay); WtS_CHANGESTAT_FN(s_nonconformity_thresholds);

WtD_CHANGESTAT_FN(d_nonconformity_thresholds); WtS_CHANGESTAT_FN(s_nonconformity_thresholds);

#endif
