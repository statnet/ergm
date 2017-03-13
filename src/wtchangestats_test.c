#include"wtchangestats_test.h"

WtD_CHANGESTAT_FN(d_test_abs_sum_minus_5){
  GET_STORAGE(double, stored_sum_ptr);
  double sum = *stored_sum_ptr;
  EXEC_THROUGH_TOGGLES({
    CHANGE_STAT[0] -= fabs(sum-5);
    sum += NEWWT-OLDWT;
    CHANGE_STAT[0] += fabs(sum-5);
    });
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
  *sum = weight-GETWT(tail, head);
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

WtD_CHANGESTAT_FN(d_test_abs_sum_minus_5_no_s){d_test_abs_sum_minus_5(ntoggles, tails, heads, weights, mtp, nwp);}
WtI_CHANGESTAT_FN(i_test_abs_sum_minus_5_no_s){i_test_abs_sum_minus_5(mtp, nwp);}
WtU_CHANGESTAT_FN(u_test_abs_sum_minus_5_no_s){u_test_abs_sum_minus_5(tail, head, weight, mtp, nwp);}



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


WtD_CHANGESTAT_FN(d_dsociomatrix){
  GET_AUX_STORAGE(double *, sm);
  
  ZERO_ALL_CHANGESTATS();
  EXEC_THROUGH_TOGGLES({
      Dyad pos = TAIL-1 + (HEAD-1)*N_NODES;
      CHANGE_STAT[pos] += NEWWT - sm[TAIL][HEAD];
    });  
}
