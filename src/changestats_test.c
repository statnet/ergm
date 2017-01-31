#include"changestats_test.h"

D_CHANGESTAT_FN(d_test_abs_edges_minus_5){
  GET_STORAGE(Edge, stored_edges_ptr);
  long int edges = *stored_edges_ptr;
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

I_CHANGESTAT_FN(i_test_abs_edges_minus_5){
  ALLOC_STORAGE(1, Edge, edges);

  *edges = N_EDGES; // Pretend this takes a long time.
}

// Finalizer is not needed.

U_CHANGESTAT_FN(u_test_abs_edges_minus_5){
  GET_STORAGE(Edge, edges);
  int i;
  FOR_EACH_TOGGLE(i){
    unsigned int edgeflag = IS_OUTEDGE(TAIL(i), HEAD(i));
    *edges += edgeflag ? - 1 : 1;
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

S_CHANGESTAT_FN(s_test_abs_edges_minus_5){
  GET_STORAGE(Edge, edges);
  // Storage uninitialized: compute from scratch.
  if(!edges){
    CHANGE_STAT[0] = labs((long int)N_EDGES-5);
  }else{ // Storage initialized: use it.
    CHANGE_STAT[0] = labs((long int)*edges-5);
  }
}

D_CHANGESTAT_FN(d_test_abs_edges_minus_5_no_s){d_test_abs_edges_minus_5(ntoggles, tails, heads, mtp, nwp);}
I_CHANGESTAT_FN(i_test_abs_edges_minus_5_no_s){i_test_abs_edges_minus_5(mtp, nwp);}
U_CHANGESTAT_FN(u_test_abs_edges_minus_5_no_s){u_test_abs_edges_minus_5(ntoggles, tails, heads, mtp, nwp);}

I_CHANGESTAT_FN(i__isociomatrix){
  // Note: the following code first sets up a 2D array indexed from 0, then shifts all pointers by -1 so that sm[t][h] would work for vertex IDs. 
  
  ALLOC_AUX_STORAGE(N_TAILS, int*, sm); // An int ** vector of pointers to ints in slot owned by this auxiliary.
  
  Dyad sm_size = BIPARTITE? N_TAILS*N_HEADS : DIRECTED ? N_NODES*N_NODES : N_NODES*(N_NODES+1)/2; // For consistency, and possible future capabilities, include diagonal.
  ALLOC_STORAGE(sm_size, int, data); // An int * to data.

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
  for(Vertex tail=1; tail <= N_TAILS; tail++){
    Vertex head;
    Edge e;
    STEP_THROUGH_OUTEDGES(tail, e, head) {
      sm[tail][head] = 1;
    }
  }
}

U_CHANGESTAT_FN(u__isociomatrix){
  GET_AUX_STORAGE(int*, sm);
  int i;
  FOR_EACH_TOGGLE(i){
    Vertex tail=TAIL(i), head=HEAD(i);
    sm[tail][head]  = 1 - sm[tail][head];
  }
}

F_CHANGESTAT_FN(f__isociomatrix){
  // If we hadn't shifted the pointers by -1, this would not have been necessary.
  GET_AUX_STORAGE(double*, sm);
  free(sm + 1);
  nwp->aux_storage[(unsigned int) INPUT_PARAM[0]] = NULL;
  // nwp->storage was not shifted, so it'll be freed automatically.
}


D_CHANGESTAT_FN(d_isociomatrix){
  GET_AUX_STORAGE(int *, sm);
  int i;
  
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i){
    Vertex tail=TAIL(i), head=HEAD(i);
    Dyad pos = tail-1 + (head-1)*N_NODES;
    int edgeflag = sm[tail][head];
    if(CHANGE_STAT[pos]==0)
      CHANGE_STAT[pos] += edgeflag? -1 : +1;
    else CHANGE_STAT[pos] *= -1;
  }
}
