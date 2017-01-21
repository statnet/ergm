#ifndef CHANGESTATS_TEST_H
#define CHANGESTATS_TEST_H

#include "edgetree.h"
#include "changestat.h"

D_CHANGESTAT_FN(d_test_abs_edges_minus_5);
U_CHANGESTAT_FN(u_test_abs_edges_minus_5);
S_CHANGESTAT_FN(s_test_abs_edges_minus_5);

D_CHANGESTAT_FN(d_test_abs_edges_minus_5_no_s);
U_CHANGESTAT_FN(u_test_abs_edges_minus_5_no_s);

#define ISOCIOMATRIX_CELL(sm, tail, head) (*(((int *) sm) + tail-1 + (head-1)*N_NODES))

U_CHANGESTAT_FN(u__isociomatrix);
D_CHANGESTAT_FN(d_isociomatrix);

#endif
