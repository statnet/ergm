/*  File src/wtchangestats_test.c in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2022 Statnet Commons
 */
#include "ergm_wtedgetree.h"
#include "ergm_wtchangestat.h"
#include "ergm_storage.h"

WtC_CHANGESTAT_FN(c_test_abs_sum_minus_5){
  GET_STORAGE(double, stored_sum_ptr);
  double sum = *stored_sum_ptr;
    CHANGE_STAT[0] = -fabs(sum-5);
    CHANGE_STAT[0] += fabs(sum-5 + weight - edgestate);
}

WtI_CHANGESTAT_FN(i_test_abs_sum_minus_5){
  ALLOC_STORAGE(1, double, sum);
  *sum = 0;
  EXEC_THROUGH_NET_EDGES(tail, e1, head, y, {
      *sum+=y;
      head=head; e1=e1; // Prevent a compiler warning.
    });
}

WtU_CHANGESTAT_FN(u_test_abs_sum_minus_5){
  GET_STORAGE(double, sum);
  *sum += weight-edgestate;
}

WtS_CHANGESTAT_FN(s_test_abs_sum_minus_5){
  GET_STORAGE(double, sum);
  // Storage uninitialized: compute from scratch.
  if(!sum){
    double sum = 0;
    EXEC_THROUGH_NET_EDGES(tail, e1, head, y, {
	sum+=y;
	head=head; e1=e1; // Prevent a compiler warning.
      });
    CHANGE_STAT[0] = fabs(sum-5);
  }else{ // Storage initialized: use it.
    CHANGE_STAT[0] = fabs(*sum-5);
  }
}

WtC_CHANGESTAT_FN(c_test_abs_sum_minus_5_no_s){c_test_abs_sum_minus_5(tail, head, weight, mtp, nwp, edgestate);}
WtI_CHANGESTAT_FN(i_test_abs_sum_minus_5_no_s){i_test_abs_sum_minus_5(mtp, nwp);}
WtU_CHANGESTAT_FN(u_test_abs_sum_minus_5_no_s){u_test_abs_sum_minus_5(tail, head, weight, mtp, nwp, edgestate);}



WtI_CHANGESTAT_FN(i__dsociomatrix){
  ALLOC_AUX_SOCIOMATRIX(double, sm);
  
  // Now, populate the sociomatrix.
  EXEC_THROUGH_NET_EDGES(t, h, e, w, {
      sm[t][h] = w;
    });
}

WtU_CHANGESTAT_FN(u__dsociomatrix){
  GET_AUX_STORAGE(double*, sm);
  sm[tail][head] = weight;
}

WtF_CHANGESTAT_FN(f__dsociomatrix){
  FREE_AUX_SOCIOMATRIX;
}


WtC_CHANGESTAT_FN(c_dsociomatrix){
  GET_AUX_STORAGE(double *, sm);
  
      Dyad pos = tail-1 + (head-1)*N_NODES;
      CHANGE_STAT[pos] = weight - sm[tail][head];
}


WtI_CHANGESTAT_FN(i__sum){
  ALLOC_AUX_STORAGE(1, double, sum);
  *sum = 0;
  EXEC_THROUGH_NET_EDGES(tail, e1, head, y, {
      *sum+=y;
      (void) head; (void) e1; // Prevent a compiler warning.
    });
}

WtU_CHANGESTAT_FN(u__sum){
  GET_AUX_STORAGE(double, sum);
  *sum += weight-edgestate;
}

WtC_CHANGESTAT_FN(c_test_abs_sum_minus_5_aux){
  GET_AUX_STORAGE(double, stored_sum_ptr);
  double sum = *stored_sum_ptr;
    CHANGE_STAT[0] = -fabs(sum-5);
    CHANGE_STAT[0] += fabs(sum-5 + weight - edgestate);
}
