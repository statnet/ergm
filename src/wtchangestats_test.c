#include"wtchangestats_test.h"

WtD_CHANGESTAT_FN(d_test_abs_sum_minus_5){
  double sum = *(double *)mtp->storage;
  ZERO_ALL_CHANGESTATS();
  EXEC_THROUGH_TOGGLES({
    CHANGE_STAT[0] -= fabs(sum-5);
    sum += NEWWT-OLDWT;
    CHANGE_STAT[0] += fabs(sum-5);
    });
}

WtU_CHANGESTAT_FN(u_test_abs_sum_minus_5){
  double *sum;
  // Uninitialized
  if(!mtp->storage){
    // Allocate and calculate the contents of the storage space. Note
    // that it's OK to allocate and deallocate temporary storage here
    // using, say, R_alloc(), since this code segment will not run very often.
    mtp->storage = malloc(sizeof(double));
    // Pretend calculating the number of sum takes a long, long time.
    sum = (double *)mtp->storage;
    *sum = 0;
    for (Vertex tail=1; tail <= N_NODES; tail++){
      EXEC_THROUGH_FOUTEDGES(tail, e1, head, y, {
	  *sum+=y;
	  head=head; // Prevent a compiler warning.
	});
    }
  }else sum = (double *)mtp->storage; 
  // Note that we need to check if there are any toggles to be applied whether or not we just initialized.
  EXEC_THROUGH_TOGGLES({
      *sum += NEWWT-OLDWT;
    });
}

WtS_CHANGESTAT_FN(s_test_abs_sum_minus_5){
  // Storage uninitialized: compute from scratch.
  if(!mtp->storage){
    double sum = 0;
    for (Vertex tail=1; tail <= N_NODES; tail++){
      EXEC_THROUGH_FOUTEDGES(tail, e1, head, y, {
	  sum+=y;
	  head=head; // Prevent a compiler warning.
	});
    }
    CHANGE_STAT[0] = fabs((long int)sum-5);
  }else{ // Storage initialized: use it.
    CHANGE_STAT[0] = fabs(*(double *)mtp->storage-5);
  }
}

WtD_CHANGESTAT_FN(d_test_abs_sum_minus_5_no_s){d_test_abs_sum_minus_5(ntoggles, tails, heads, weights, mtp, nwp);}
WtU_CHANGESTAT_FN(u_test_abs_sum_minus_5_no_s){u_test_abs_sum_minus_5(ntoggles, tails, heads, weights, mtp, nwp);}
