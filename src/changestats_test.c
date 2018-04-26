/*  File src/changestats_test.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2017 Statnet Commons
 */
#include"changestats_test.h"
#include"ergm_changestats_auxnet.h"

C_CHANGESTAT_FN(c_test_abs_edges_minus_5){
  GET_STORAGE(Edge, stored_edges_ptr);
  long int edges = *stored_edges_ptr;
  unsigned int edgeflag = IS_OUTEDGE(tail, head);
  CHANGE_STAT[0] = -labs(edges-5);
  CHANGE_STAT[0] += labs(edges-5 + (edgeflag?-1:1));
}

I_CHANGESTAT_FN(i_test_abs_edges_minus_5){
  ALLOC_STORAGE(1, Edge, edges);

  *edges = N_EDGES; // Pretend this takes a long time.
}

// Finalizer is not needed.

U_CHANGESTAT_FN(u_test_abs_edges_minus_5){
  GET_STORAGE(Edge, edges);
  *edges += IS_OUTEDGE(tail, head) ? - 1 : 1;
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

C_CHANGESTAT_FN(c_test_abs_edges_minus_5_no_s){c_test_abs_edges_minus_5(tail, head, mtp, nwp);}
I_CHANGESTAT_FN(i_test_abs_edges_minus_5_no_s){i_test_abs_edges_minus_5(mtp, nwp);}
U_CHANGESTAT_FN(u_test_abs_edges_minus_5_no_s){u_test_abs_edges_minus_5(tail, head, mtp, nwp);}

C_CHANGESTAT_FN(c_isociomatrix){
  GET_AUX_STORAGE(int *, sm);
  
  ZERO_ALL_CHANGESTATS(i);
    Dyad pos = tail-1 + (head-1)*N_NODES;
    CHANGE_STAT[pos] += sm[tail][head]? -1 : +1;
}

C_CHANGESTAT_FN(c_disc_inter_union_net){
  GET_AUX_STORAGE_NUM(Network, dnwp, 0);
  GET_AUX_STORAGE_NUM(Network, inwp, 1);
  GET_AUX_STORAGE_NUM(Network, unwp, 2);

  int nwedge = IS_OUTEDGE(tail, head)!=0;
  int refedge = dEdgeListSearch(tail, head, INPUT_PARAM+3)!=0;
  
  CHANGE_STAT[0] = nwedge!=refedge ? -1 : +1;
  CHANGE_STAT[1] = refedge ? (nwedge ? -1 : +1) : 0;
  CHANGE_STAT[2] = !refedge ? (nwedge ? -1 : +1) : 0;

  CHANGE_STAT[3] = (dnwp->nedges+CHANGE_STAT[0])*(dnwp->nedges+CHANGE_STAT[0]) - dnwp->nedges*dnwp->nedges;
  CHANGE_STAT[4] = (inwp->nedges+CHANGE_STAT[1])*(inwp->nedges+CHANGE_STAT[1]) - inwp->nedges*inwp->nedges;
  CHANGE_STAT[5] = (unwp->nedges+CHANGE_STAT[2])*(unwp->nedges+CHANGE_STAT[2]) - unwp->nedges*unwp->nedges;
}

/*****************
 changestat: c__edges_times
*****************/
C_CHANGESTAT_FN(c__edges_times) {
  int edgeflag;
  edgeflag = IS_OUTEDGE(tail, head);
  CHANGE_STAT[0] += edgeflag ? - *INPUT_PARAM : *INPUT_PARAM;
}

S_CHANGESTAT_FN(s__edges_tests) {
  CHANGE_STAT[0] = N_EDGES * *INPUT_PARAM;
}
