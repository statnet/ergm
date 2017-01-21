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
  INIT_STORAGE(Edge, edges, {
      *edges = N_EDGES; // Pretend this takes a long time.
    });

  int i;
  FOR_EACH_TOGGLE(i){
    unsigned int edgeflag = IS_OUTEDGE(TAIL(i), HEAD(i));
    *edges += edgeflag ? - 1 : 1;
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


U_CHANGESTAT_FN(u__isociomatrix){
  unsigned int myslot = INPUT_PARAM[0];
  int *sm = nwp->aux_storage[myslot];
  if(!sm){
    sm = (int *) (nwp->aux_storage[myslot] = malloc(sizeof(int)*N_NODES*N_NODES));
    memset(sm, 0, sizeof(int)*N_NODES*N_NODES);
    for(Vertex tail=1; tail <= N_NODES; tail++){
      Vertex head;
      Edge e;
      STEP_THROUGH_OUTEDGES(tail, e, head) {
	ISOCIOMATRIX_CELL(sm, tail, head) = 1; 
      }
    }
  }

  int i;
  FOR_EACH_TOGGLE(i){
    Vertex tail=TAIL(i), head=HEAD(i);
    ISOCIOMATRIX_CELL(sm, tail, head) = 1 - ISOCIOMATRIX_CELL(sm, tail, head);
  }  
}

D_CHANGESTAT_FN(d_isociomatrix){
  unsigned int slot = INPUT_PARAM[0];
  int *sm = nwp->aux_storage[slot];
  int i;
  
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i){
    Vertex tail=TAIL(i), head=HEAD(i);
    Dyad pos = tail-1 + (head-1)*N_NODES;
    int edgeflag = ISOCIOMATRIX_CELL(sm, tail, head);
    CHANGE_STAT[pos] += edgeflag? -1 : +1;
  }  
}
