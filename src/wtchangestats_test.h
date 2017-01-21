#ifndef WTCHANGESTATS_TEST_H
#define WTCHANGESTATS_TEST_H

#include "wtedgetree.h"
#include "wtchangestat.h"

WtD_CHANGESTAT_FN(d_test_abs_sum_minus_5);
WtU_CHANGESTAT_FN(u_test_abs_sum_minus_5);
WtS_CHANGESTAT_FN(s_test_abs_sum_minus_5);

WtD_CHANGESTAT_FN(d_test_abs_sum_minus_5_no_s);
WtU_CHANGESTAT_FN(u_test_abs_sum_minus_5_no_s);

#define DSOCIOMATRIX_CELL(sm, tail, head) (*(((double *) sm) + tail-1 + (head-1)*N_NODES))

WtU_CHANGESTAT_FN(u__dsociomatrix);
WtD_CHANGESTAT_FN(d_dsociomatrix);

#endif
