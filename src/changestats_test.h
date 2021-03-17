/*  File src/changestats_test.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2017 Statnet Commons
 */
#ifndef CHANGESTATS_TEST_H
#define CHANGESTATS_TEST_H

#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_storage.h"

D_CHANGESTAT_FN(d_test_abs_edges_minus_5);
I_CHANGESTAT_FN(i_test_abs_edges_minus_5);
U_CHANGESTAT_FN(u_test_abs_edges_minus_5);
F_CHANGESTAT_FN(f_test_abs_edges_minus_5);
S_CHANGESTAT_FN(s_test_abs_edges_minus_5);

C_CHANGESTAT_FN(c_test_abs_edges_minus_5_no_s);
I_CHANGESTAT_FN(i_test_abs_edges_minus_5_no_s);
U_CHANGESTAT_FN(u_test_abs_edges_minus_5_no_s);

C_CHANGESTAT_FN(c_isociomatrix);

C_CHANGESTAT_FN(c_disc_inter_union_net);


#endif
