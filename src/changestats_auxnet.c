/*  File src/changestats_test.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2017 Statnet Commons
 */
#include "ergm_changestats_auxnet.h"
#include "ergm_dyad_hashmap.h"

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
