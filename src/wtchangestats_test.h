/*  File src/wtchangestats_test.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2017 Statnet Commons
 */
#ifndef WTCHANGESTATS_TEST_H
#define WTCHANGESTATS_TEST_H

#include "ergm_wtedgetree.h"
#include "ergm_wtchangestat.h"
#include "ergm_storage.h"

WtD_CHANGESTAT_FN(d_test_abs_sum_minus_5);
WtI_CHANGESTAT_FN(i_test_abs_sum_minus_5);
WtU_CHANGESTAT_FN(u_test_abs_sum_minus_5);
WtF_CHANGESTAT_FN(f_test_abs_sum_minus_5);
WtS_CHANGESTAT_FN(s_test_abs_sum_minus_5);

WtD_CHANGESTAT_FN(d_test_abs_sum_minus_5_no_s);
WtI_CHANGESTAT_FN(i_test_abs_sum_minus_5_no_s);
WtU_CHANGESTAT_FN(u_test_abs_sum_minus_5_no_s);

WtI_CHANGESTAT_FN(i__dsociomatrix);
WtU_CHANGESTAT_FN(u__dsociomatrix);
WtF_CHANGESTAT_FN(f__dsociomatrix);

WtD_CHANGESTAT_FN(d_dsociomatrix);

#endif
