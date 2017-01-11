#include"changestats_test.h"

D_CHANGESTAT_FN(d_test_abs_edges_minus_5){
  long int edges = *(Edge *)mtp->storage;
  int i;
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i){
    unsigned int edgeflag = IS_OUTEDGE(TAIL(i), HEAD(i));
    CHANGE_STAT[0] -= labs(edges-5);
    edges += edgeflag ? - 1 : 1;
    CHANGE_STAT[0] += labs(edges-5);
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

U_CHANGESTAT_FN(u_test_abs_edges_minus_5){
  // Uninitialized
  if(!mtp->storage){
    // Allocate and calculate the contents of the storage space. Note
    // that it's OK to allocate and deallocate temporary storage here
    // using, say, R_alloc(), since this code segment will not run very often.
    mtp->storage = malloc(sizeof(Edge));
    // Pretend calculating the number of edges takes a long, long time.
    *(Edge *)mtp->storage = N_EDGES; // This is a bit convoluted: cast the void * to an Edge *, then dereference and assign.
  }
  // Note that we need to check if there are any toggles to be applied whether or not we just initialized.

  int i;
  FOR_EACH_TOGGLE(i){
    unsigned int edgeflag = IS_OUTEDGE(TAIL(i), HEAD(i));
    *((Edge *)mtp->storage) += edgeflag ? - 1 : 1;
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

S_CHANGESTAT_FN(s_test_abs_edges_minus_5){
  // Storage uninitialized: compute from scratch.
  if(!mtp->storage){
    CHANGE_STAT[0] = labs((long int)N_EDGES-5);
  }else{ // Storage initialized: use it.
    CHANGE_STAT[0] = labs((long int)*(Edge *)mtp->storage-5);
  }
}

D_CHANGESTAT_FN(d_test_abs_edges_minus_5_no_s){d_test_abs_edges_minus_5(ntoggles, tails, heads, mtp, nwp);}
U_CHANGESTAT_FN(u_test_abs_edges_minus_5_no_s){u_test_abs_edges_minus_5(ntoggles, tails, heads, mtp, nwp);}
