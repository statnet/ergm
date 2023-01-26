/*  File src/changestats_auxnetalike.c in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2023 Statnet Commons
 */
#include "ergm_changestats_auxnet.h"
#include "ergm_dyad_hashmap.h"
#include "ergm_dyad_hashmap_utils.h"

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

I_CHANGESTAT_FN(i__discord_isociomatrix){
  ALLOC_AUX_SOCIOMATRIX(int, sm);
  GET_AUX_STORAGE_NUM(StoreAuxnet, storage, 1);

  nwp = storage->onwp; // So that we can use the macros.
  // Now, populate the sociomatrix.
  EXEC_THROUGH_NET_EDGES(tail, head, e, {
      sm[tail][head] = 1;
    });
}

U_CHANGESTAT_FN(u__discord_isociomatrix){
  GET_AUX_STORAGE(int*, sm);
  sm[tail][head]  = 1 - sm[tail][head];
}

F_CHANGESTAT_FN(f__discord_isociomatrix){
  FREE_AUX_SOCIOMATRIX;
}

I_CHANGESTAT_FN(i__discord_net_DyadSet){
  ALLOC_AUX_STORAGE(1, StoreDyadSetAndRefEL, storage);
  StoreDyadSet *dnwp = storage->nwp = NetworkToDyadSet(nwp);
  int *ref_el = storage->ref_el = IINPUT_PARAM;
  
  Edge nedges = *ref_el;
  for(Edge i=0; i<nedges; i++){
    Vertex tail=ref_el[1+i], head=ref_el[1+nedges+i];
    DyadSetToggle(tail,head, dnwp);
  }
}

U_CHANGESTAT_FN(u__discord_net_DyadSet){
  GET_AUX_STORAGE(StoreDyadSetAndRefEL, storage);
  StoreDyadSet *dnwp = storage->nwp;

  DyadSetToggle(tail,head, dnwp);
}

F_CHANGESTAT_FN(f__discord_net_DyadSet){
  GET_AUX_STORAGE(StoreDyadSetAndRefEL, storage);
  StoreDyadSet *dnwp = storage->nwp;

  kh_destroy(DyadSet, dnwp);
}

I_CHANGESTAT_FN(i__intersect_net_DyadSet){
  ALLOC_AUX_STORAGE(1, StoreDyadSetAndRefEL, storage);
  StoreDyadSet *dnwp = storage->nwp = kh_init(DyadSet); dnwp->directed=DIRECTED;
  int *ref_el = storage->ref_el = IINPUT_PARAM;
  
  Edge nedges = *ref_el;
  for(Edge i=0; i<nedges; i++){
    Vertex tail=ref_el[1+i], head=ref_el[1+nedges+i];
    if(IS_OUTEDGE(tail, head)) {
      DyadSetToggle(tail,head, dnwp);
    }
  }
}

U_CHANGESTAT_FN(u__intersect_net_DyadSet){
  GET_AUX_STORAGE(StoreDyadSetAndRefEL, storage);
  StoreDyadSet *dnwp = storage->nwp;
  int *ref_el = storage->ref_el;
  // only toggle if the edge is in y0. otherwise changing y1 won't matter.
  if(iEdgeListSearch(tail, head, ref_el))
    DyadSetToggle(tail,head, dnwp);
}

F_CHANGESTAT_FN(f__intersect_net_DyadSet){
  GET_AUX_STORAGE(StoreDyadSetAndRefEL, storage);
  StoreDyadSet *dnwp = storage->nwp;

  kh_destroy(DyadSet, dnwp);
}

I_CHANGESTAT_FN(i__intersect_net_toggles_in_list_DyadSet){
  ALLOC_AUX_STORAGE(1, StoreDyadSetAndRefEL, storage);
  StoreDyadSet *dnwp = storage->nwp = kh_init(DyadSet); dnwp->directed=DIRECTED;
  int *ref_el = storage->ref_el = IINPUT_PARAM;
  
  Edge nedges = *ref_el;
  for(Edge i=0; i<nedges; i++){
    Vertex tail=ref_el[1+i], head=ref_el[1+nedges+i];
    if(IS_OUTEDGE(tail, head)!=0)
      DyadSetToggle(tail,head, dnwp);
  }
}

U_CHANGESTAT_FN(u__intersect_net_toggles_in_list_DyadSet){
  GET_AUX_STORAGE(StoreDyadSetAndRefEL, storage);
  StoreDyadSet *dnwp = storage->nwp;

  DyadSetToggle(tail,head, dnwp);
}

F_CHANGESTAT_FN(f__intersect_net_toggles_in_list_DyadSet){
  GET_AUX_STORAGE(StoreDyadSetAndRefEL, storage);
  StoreDyadSet *dnwp = storage->nwp;

  kh_destroy(DyadSet, dnwp);
}

I_CHANGESTAT_FN(i__union_net_DyadSet){
  ALLOC_AUX_STORAGE(1, StoreDyadSetAndRefEL, storage);
  StoreDyadSet *dnwp = storage->nwp = NetworkToDyadSet(nwp);
  int *ref_el = storage->ref_el = IINPUT_PARAM;
  
  Edge nedges = *ref_el;
  for(Edge i=0; i<nedges; i++){
    Vertex tail=ref_el[1+i], head=ref_el[1+nedges+i];
    if(IS_OUTEDGE(tail, head)==0)
      DyadSetToggle(tail,head, dnwp);
  }
}

U_CHANGESTAT_FN(u__union_net_DyadSet){
  GET_AUX_STORAGE(StoreDyadSetAndRefEL, storage);
  StoreDyadSet *dnwp = storage->nwp;
  int *ref_el = storage->ref_el;
  // If the edge is in y0, changing y1 won't matter.
  if(iEdgeListSearch(tail, head, ref_el)==0)
    DyadSetToggle(tail,head, dnwp);
}

F_CHANGESTAT_FN(f__union_net_DyadSet){
  GET_AUX_STORAGE(StoreDyadSetAndRefEL, storage);
  StoreDyadSet *dnwp = storage->nwp;

  kh_destroy(DyadSet, dnwp);
}
