/*  File src/changestats_test.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2017 Statnet Commons
 */
#include "changestats_test.h"
#include "ergm_dyad_hashmap.h"

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

I_CHANGESTAT_FN(i__isociomatrix){
  ALLOC_AUX_SOCIOMATRIX(int, sm);
  
  // Now, populate the sociomatrix.
  EXEC_THROUGH_NET_EDGES(tail, head, e, {
      sm[tail][head] = 1;
    });
}

U_CHANGESTAT_FN(u__isociomatrix){
  GET_AUX_STORAGE(int*, sm);
  sm[tail][head]  = 1 - sm[tail][head];
}

F_CHANGESTAT_FN(f__isociomatrix){
  FREE_AUX_SOCIOMATRIX;
}


C_CHANGESTAT_FN(c_isociomatrix){
  GET_AUX_STORAGE(int *, sm);
  
  ZERO_ALL_CHANGESTATS(i);
    Dyad pos = tail-1 + (head-1)*N_NODES;
    CHANGE_STAT[pos] += sm[tail][head]? -1 : +1;
}

I_CHANGESTAT_FN(i__discord_net_Network){
  Network *dnwp = AUX_STORAGE = NetworkCopy(nwp);
  Edge nedges = INPUT_PARAM[1];
  for(Edge i=0; i<nedges; i++){
    Vertex tail=INPUT_PARAM[2+i], head=INPUT_PARAM[2+nedges+i];
    ToggleEdge(tail,head, dnwp);
  }
}

U_CHANGESTAT_FN(u__discord_net_Network){
  Network *dnwp = AUX_STORAGE;

  ToggleEdge(tail, head, dnwp);
}

F_CHANGESTAT_FN(f__discord_net_Network){
  Network *dnwp = AUX_STORAGE;

  NetworkDestroy(dnwp);
  AUX_STORAGE = NULL;
}

I_CHANGESTAT_FN(i__intersect_net_Network){
  Network *dnwp = AUX_STORAGE = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, 0, 0, NULL);
  Edge nedges = INPUT_PARAM[1];
  for(Edge i=0; i<nedges; i++){
    Vertex tail=INPUT_PARAM[2+i], head=INPUT_PARAM[2+nedges+i];
    if(IS_OUTEDGE(tail, head)!=0)
      ToggleEdge(tail,head, dnwp);
  }
}

U_CHANGESTAT_FN(u__intersect_net_Network){
  Network *dnwp = AUX_STORAGE;

  if(dEdgeListSearch(tail, head, INPUT_PARAM+1))
    ToggleEdge(tail, head, dnwp);
}

F_CHANGESTAT_FN(f__intersect_net_Network){
  Network *dnwp = AUX_STORAGE;

  NetworkDestroy(dnwp);
  AUX_STORAGE = NULL;
}

I_CHANGESTAT_FN(i__intersect_net_toggles_in_list_Network){
  Network *dnwp = AUX_STORAGE = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, 0, 0, NULL);
  Edge nedges = INPUT_PARAM[1];
  for(Edge i=0; i<nedges; i++){
    Vertex tail=INPUT_PARAM[2+i], head=INPUT_PARAM[2+nedges+i];
    if(IS_OUTEDGE(tail, head)!=0)
      ToggleEdge(tail,head, dnwp);
  }
}

U_CHANGESTAT_FN(u__intersect_net_toggles_in_list_Network){
  Network *dnwp = AUX_STORAGE;

  ToggleEdge(tail, head, dnwp);
}

F_CHANGESTAT_FN(f__intersect_net_toggles_in_list_Network){
  Network *dnwp = AUX_STORAGE;

  NetworkDestroy(dnwp);
  AUX_STORAGE = NULL;
}

I_CHANGESTAT_FN(i__union_net_Network){
  Network *dnwp = AUX_STORAGE = NetworkCopy(nwp);
  Edge nedges = INPUT_PARAM[1];
  for(Edge i=0; i<nedges; i++){
    Vertex tail=INPUT_PARAM[2+i], head=INPUT_PARAM[2+nedges+i];
    if(IS_OUTEDGE(tail, head)==0)
      ToggleEdge(tail,head, dnwp);
  }
}

U_CHANGESTAT_FN(u__union_net_Network){
  Network *dnwp = AUX_STORAGE;

  if(dEdgeListSearch(tail, head, INPUT_PARAM+1)==0)
    ToggleEdge(tail, head, dnwp);
}

F_CHANGESTAT_FN(f__union_net_Network){
  Network *dnwp = AUX_STORAGE;

  NetworkDestroy(dnwp);
  AUX_STORAGE = NULL;
}

C_CHANGESTAT_FN(c_disc_inter_union_net_Network){
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

I_CHANGESTAT_FN(i__discord_net_DyadSet){
  StoreDyadSet *dnwp = AUX_STORAGE = NetworkToDyadSet(nwp);
  Edge nedges = INPUT_PARAM[1];
  for(Edge i=0; i<nedges; i++){
    Vertex tail=INPUT_PARAM[2+i], head=INPUT_PARAM[2+nedges+i];
    DyadSetToggle(TH(tail,head,DIRECTED), dnwp);
  }
}

U_CHANGESTAT_FN(u__discord_net_DyadSet){
  StoreDyadSet *dnwp = AUX_STORAGE;

  DyadSetToggle(TH(tail,head,DIRECTED), dnwp);
}

F_CHANGESTAT_FN(f__discord_net_DyadSet){
  StoreDyadSet *dnwp = AUX_STORAGE;

  kh_destroy(DyadSet, dnwp);
  AUX_STORAGE = NULL;
}

I_CHANGESTAT_FN(i__intersect_net_DyadSet){
  StoreDyadSet *dnwp = AUX_STORAGE = kh_init(DyadSet);
  Edge nedges = INPUT_PARAM[1];
  for(Edge i=0; i<nedges; i++){
    Vertex tail=INPUT_PARAM[2+i], head=INPUT_PARAM[2+nedges+i];
    if(IS_OUTEDGE(tail, head)!=0)
      DyadSetToggle(TH(tail,head,DIRECTED), dnwp);
  }
}

U_CHANGESTAT_FN(u__intersect_net_DyadSet){
  StoreDyadSet *dnwp = AUX_STORAGE;

  if(dEdgeListSearch(tail, head, INPUT_PARAM+1))
    DyadSetToggle(TH(tail,head,DIRECTED), dnwp);
}

F_CHANGESTAT_FN(f__intersect_net_DyadSet){
  StoreDyadSet *dnwp = AUX_STORAGE;

  kh_destroy(DyadSet, dnwp);
  AUX_STORAGE = NULL;
}

I_CHANGESTAT_FN(i__intersect_net_toggles_in_list_DyadSet){
  StoreDyadSet *dnwp = AUX_STORAGE = kh_init(DyadSet);
  Edge nedges = INPUT_PARAM[1];
  for(Edge i=0; i<nedges; i++){
    Vertex tail=INPUT_PARAM[2+i], head=INPUT_PARAM[2+nedges+i];
    if(IS_OUTEDGE(tail, head)!=0)
      DyadSetToggle(TH(tail,head,DIRECTED), dnwp);
  }
}

U_CHANGESTAT_FN(u__intersect_net_toggles_in_list_DyadSet){
  StoreDyadSet *dnwp = AUX_STORAGE;

  DyadSetToggle(TH(tail,head,DIRECTED), dnwp);
}

F_CHANGESTAT_FN(f__intersect_net_toggles_in_list_DyadSet){
  StoreDyadSet *dnwp = AUX_STORAGE;

  kh_destroy(DyadSet, dnwp);
  AUX_STORAGE = NULL;
}

I_CHANGESTAT_FN(i__union_net_DyadSet){
  StoreDyadSet *dnwp = AUX_STORAGE = NetworkToDyadSet(nwp);
  Edge nedges = INPUT_PARAM[1];
  for(Edge i=0; i<nedges; i++){
    Vertex tail=INPUT_PARAM[2+i], head=INPUT_PARAM[2+nedges+i];
    if(IS_OUTEDGE(tail, head)==0)
      DyadSetToggle(TH(tail,head,DIRECTED), dnwp);
  }
}

U_CHANGESTAT_FN(u__union_net_DyadSet){
  StoreDyadSet *dnwp = AUX_STORAGE;

  if(dEdgeListSearch(tail, head, INPUT_PARAM+1)==0)
    DyadSetToggle(TH(tail,head,DIRECTED), dnwp);
}

F_CHANGESTAT_FN(f__union_net_DyadSet){
  StoreDyadSet *dnwp = AUX_STORAGE;

  kh_destroy(DyadSet, dnwp);
  AUX_STORAGE = NULL;
}

C_CHANGESTAT_FN(c_disc_inter_union_net_DyadSet){
  GET_AUX_STORAGE_NUM(StoreDyadSet, dnwp, 0);
  GET_AUX_STORAGE_NUM(StoreDyadSet, inwp, 1);
  GET_AUX_STORAGE_NUM(StoreDyadSet, unwp, 2);

  int nwedge = IS_OUTEDGE(tail, head)!=0;
  int refedge = dEdgeListSearch(tail, head, INPUT_PARAM+3)!=0;
  
  CHANGE_STAT[0] = nwedge!=refedge ? -1 : +1;
  CHANGE_STAT[1] = refedge ? (nwedge ? -1 : +1) : 0;
  CHANGE_STAT[2] = !refedge ? (nwedge ? -1 : +1) : 0;

  CHANGE_STAT[3] = (dnwp->size+CHANGE_STAT[0])*(dnwp->size+CHANGE_STAT[0]) - dnwp->size*dnwp->size;
  CHANGE_STAT[4] = (inwp->size+CHANGE_STAT[1])*(inwp->size+CHANGE_STAT[1]) - inwp->size*inwp->size;
  CHANGE_STAT[5] = (unwp->size+CHANGE_STAT[2])*(unwp->size+CHANGE_STAT[2]) - unwp->size*unwp->size;
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
