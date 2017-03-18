#ifndef CHANGESTATS_TEST_H
#define CHANGESTATS_TEST_H

#include "edgetree.h"
#include "changestat.h"

D_CHANGESTAT_FN(d_test_abs_edges_minus_5);
I_CHANGESTAT_FN(i_test_abs_edges_minus_5);
U_CHANGESTAT_FN(u_test_abs_edges_minus_5);
F_CHANGESTAT_FN(f_test_abs_edges_minus_5);
S_CHANGESTAT_FN(s_test_abs_edges_minus_5);

C_CHANGESTAT_FN(c_test_abs_edges_minus_5_no_s);
I_CHANGESTAT_FN(i_test_abs_edges_minus_5_no_s);
U_CHANGESTAT_FN(u_test_abs_edges_minus_5_no_s);

I_CHANGESTAT_FN(i__isociomatrix);
U_CHANGESTAT_FN(u__isociomatrix);
F_CHANGESTAT_FN(f__isociomatrix);

C_CHANGESTAT_FN(c_isociomatrix);

I_CHANGESTAT_FN(i__discord_net);
U_CHANGESTAT_FN(u__discord_net);
F_CHANGESTAT_FN(f__discord_net);

I_CHANGESTAT_FN(i__intersect_net);
U_CHANGESTAT_FN(u__intersect_net);
F_CHANGESTAT_FN(f__intersect_net);

I_CHANGESTAT_FN(i__intersect_net_toggles_in_list);
U_CHANGESTAT_FN(u__intersect_net_toggles_in_list);
F_CHANGESTAT_FN(f__intersect_net_toggles_in_list);

I_CHANGESTAT_FN(i__union_net);
U_CHANGESTAT_FN(u__union_net);
F_CHANGESTAT_FN(f__union_net);

C_CHANGESTAT_FN(c_disc_inter_union_net);


#endif
