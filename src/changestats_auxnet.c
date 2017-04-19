#include "ergm_changestats_auxnet.h"

I_CHANGESTAT_FN(i__isociomatrix){
  ALLOC_AUX_SOCIOMATRIX(int, sm);
  
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
  sm[tail][head]  = 1 - sm[tail][head];
}

F_CHANGESTAT_FN(f__isociomatrix){
  FREE_AUX_SOCIOMATRIX;
}


I_CHANGESTAT_FN(i__discord_net){
  ALLOC_AUX_STORAGE(1, Network, dnwp);
  NetworkCopy(dnwp, nwp);
  Edge nedges = INPUT_PARAM[1];
  for(Edge i=0; i<nedges; i++){
    Vertex tail=INPUT_PARAM[2+i], head=INPUT_PARAM[2+nedges+i];
    ToggleEdge(tail,head, dnwp);
  }
}

U_CHANGESTAT_FN(u__discord_net){
  GET_AUX_STORAGE(Network, dnwp);

  ToggleEdge(tail, head, dnwp);
}

F_CHANGESTAT_FN(f__discord_net){
  GET_AUX_STORAGE(Network, dnwp);

  NetworkDestroy(dnwp);
}

I_CHANGESTAT_FN(i__intersect_net){
  ALLOC_AUX_STORAGE(1, Network, dnwp);
  *dnwp = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, 0, 0, NULL);
  Edge nedges = INPUT_PARAM[1];
  for(Edge i=0; i<nedges; i++){
    Vertex tail=INPUT_PARAM[2+i], head=INPUT_PARAM[2+nedges+i];
    if(IS_OUTEDGE(tail, head)!=0)
      ToggleEdge(tail,head, dnwp);
  }
}

U_CHANGESTAT_FN(u__intersect_net){
  GET_AUX_STORAGE(Network, dnwp);

  if(dEdgeListSearch(tail, head, INPUT_PARAM+1))
    ToggleEdge(tail, head, dnwp);
}

F_CHANGESTAT_FN(f__intersect_net){
  GET_AUX_STORAGE(Network, dnwp);

  NetworkDestroy(dnwp);
}

I_CHANGESTAT_FN(i__intersect_net_toggles_in_list){
  ALLOC_AUX_STORAGE(1, Network, dnwp);
  *dnwp = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, 0, 0, NULL);
  Edge nedges = INPUT_PARAM[1];
  for(Edge i=0; i<nedges; i++){
    Vertex tail=INPUT_PARAM[2+i], head=INPUT_PARAM[2+nedges+i];
    if(IS_OUTEDGE(tail, head)!=0)
      ToggleEdge(tail,head, dnwp);
  }
}

U_CHANGESTAT_FN(u__intersect_net_toggles_in_list){
  GET_AUX_STORAGE(Network, dnwp);

  ToggleEdge(tail, head, dnwp);
}

F_CHANGESTAT_FN(f__intersect_net_toggles_in_list){
  GET_AUX_STORAGE(Network, dnwp);

  NetworkDestroy(dnwp);
}

I_CHANGESTAT_FN(i__union_net){
  ALLOC_AUX_STORAGE(1, StoreNetAndRefEL, storage);
  Network *dnwp = &(storage->nw);
  double *ref_el = storage->ref_el = INPUT_PARAM + 1;
  NetworkCopy(dnwp, nwp);
  Edge nedges = *ref_el;
  for(Edge i=0; i<nedges; i++){
    Vertex tail=ref_el[1+i], head=ref_el[1+nedges+i];
    if(IS_OUTEDGE(tail, head)==0)
      ToggleEdge(tail,head, dnwp);
  }
}

U_CHANGESTAT_FN(u__union_net){
  GET_AUX_STORAGE(StoreNetAndRefEL, storage);
  Network *dnwp = &(storage->nw);
  double *ref_el = storage->ref_el;

  if(dEdgeListSearch(tail, head, ref_el)==0)
    ToggleEdge(tail, head, dnwp);
}

F_CHANGESTAT_FN(f__union_net){
  GET_AUX_STORAGE(StoreNetAndRefEL, storage);
  Network *dnwp = &(storage->nw);

  NetworkDestroy(dnwp);
}

/* blockdiag_net:

   maintains a network that mirrors the main network but excludes ties
   not within specified blocks; also exports a numeric vector of block
   memberships
 */

I_CHANGESTAT_FN(i__blockdiag_net){
  ALLOC_AUX_STORAGE(1, StoreNetAndBID, blkinfo);
  blkinfo->nw = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, 0, 0, NULL);
  Network *bnwp = &(blkinfo->nw);
  double *b = blkinfo->b = INPUT_PARAM;

  for(Vertex tail=1; tail <= N_TAILS; tail++){
    Vertex head;
    Edge e;
    STEP_THROUGH_OUTEDGES(tail, e, head) {
      if(b[tail]==b[head])
	ToggleEdge(tail, head, bnwp);
    }
  }
}

U_CHANGESTAT_FN(u__blockdiag_net){
  GET_AUX_STORAGE(StoreNetAndBID, blkinfo);
  Network *bnwp = &(blkinfo->nw);
  double *b = blkinfo->b;

  if(b[tail]==b[head])
    ToggleEdge(tail, head, bnwp);
}

F_CHANGESTAT_FN(f__blockdiag_net){
  GET_AUX_STORAGE(StoreNetAndBID, blkinfo);
  Network *bnwp = &(blkinfo->nw);
  NetworkDestroy(bnwp);
}
