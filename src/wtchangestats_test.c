#include"wtchangestats_test.h"

WtD_CHANGESTAT_FN(d_test_abs_sum_minus_5){
  GET_STORAGE(double, stored_sum_ptr);
  double sum = *stored_sum_ptr;
  ZERO_ALL_CHANGESTATS();
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
  EXEC_THROUGH_TOGGLES({
      *sum += NEWWT-OLDWT;
    });
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
WtU_CHANGESTAT_FN(u_test_abs_sum_minus_5_no_s){u_test_abs_sum_minus_5(ntoggles, tails, heads, weights, mtp, nwp);}



WtI_CHANGESTAT_FN(i__dsociomatrix){
  // Note: the following code first sets up a 2D array indexed from 0, then shifts all pointers by -1 so that sm[t][h] would work for vertex IDs. 
  
  ALLOC_AUX_STORAGE(N_TAILS, double*, sm); // An double ** vector of pointers to ints in slot owned by this auxiliary.
  
  Dyad sm_size = BIPARTITE? N_TAILS*N_HEADS : DIRECTED ? N_NODES*N_NODES : N_NODES*(N_NODES+1)/2; // For consistency, and possible future capabilities, include diagonal.
  ALLOC_STORAGE(sm_size, double, data); // An double * to data.

  Dyad pos = 0; // Start of the next row's data in the data vector.
  for(Vertex t=0; t<N_TAILS; t++){
    // First set up the pointer to the right location in the data vector,    
    if(BIPARTITE){
      sm[t] = data+pos - N_TAILS; // This is so that sm[t][h=BIPARTITE] would address the 0th element of that row.
      pos += N_HEADS;
    }else if(DIRECTED){
      sm[t] = data+pos;
      pos += N_HEADS;
    }else{ // Undirected.
      sm[t] = data+pos - t; // tail <= head, so this is so that sm[t][h=t] would address the 0th element of that row. 
      pos += N_HEADS-t+1; // Each row has N_NODES - t + 1 elements (including diagonal).
    }
    sm[t]--; // Now, shift the pointer by -1. 
  }

  sm--; // Shift the pointer array by -1.
  nwp->aux_storage[(unsigned int) INPUT_PARAM[0]] = sm; // This is needed to make sure the pointer array itself is updated.

  // Now, populate the sociomatrix.
  EXEC_THROUGH_NET_EDGES(t, h, e, w, {
      sm[t][h] = w;
    });
}


WtU_CHANGESTAT_FN(u__dsociomatrix){
  GET_AUX_STORAGE(double*, sm);
  EXEC_THROUGH_TOGGLES({
      sm[TAIL][HEAD]  = NEWWT;
    });
}

WtF_CHANGESTAT_FN(f__dsociomatrix){
  unsigned int myslot = (unsigned int) INPUT_PARAM[0];
  // If we hadn't shifted the pointers by -1, this would not have been necessary.
  GET_AUX_STORAGE(double*, sm);
  free(sm + 1);
  nwp->aux_storage[myslot] = NULL;
  // nwp->storage was not shifted, so it'll be freed automatically.
}


WtD_CHANGESTAT_FN(d_dsociomatrix){
  GET_AUX_STORAGE(double *, sm);
  
  ZERO_ALL_CHANGESTATS();
  EXEC_THROUGH_TOGGLES({
      Dyad pos = TAIL-1 + (HEAD-1)*N_NODES;
      CHANGE_STAT[pos] += NEWWT - sm[TAIL][HEAD];
    });  
}
