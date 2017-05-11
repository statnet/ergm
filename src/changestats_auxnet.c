#include "ergm_changestats_auxnet.h"

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

// sets aux network to y0 XOR y1
I_CHANGESTAT_FN(i__discord_net){
  ALLOC_AUX_STORAGE(1, StoreNetAndRefEL, storage);
  Network *dnwp = &(storage->nw);
  double *ref_el = storage->ref_el = INPUT_PARAM + 1;
  NetworkCopy(dnwp, nwp);
  
  Edge nedges = *ref_el;
  for(Edge i=0; i<nedges; i++){
    Vertex tail=ref_el[1+i], head=ref_el[1+nedges+i];
    ToggleEdge(tail,head, dnwp);
  }
}

U_CHANGESTAT_FN(u__discord_net){
  GET_AUX_STORAGE(StoreNetAndRefEL, storage);
  Network *dnwp = &(storage->nw);
  
  ToggleEdge(tail, head, dnwp);
}

F_CHANGESTAT_FN(f__discord_net){
  GET_AUX_STORAGE(Network, dnwp);
  
  NetworkDestroy(dnwp);
}

// Initialize empty aux network. Then loop through the edges of y0 (ref_el).
// If the edge also exists in y1, then toggle it off in dnwp.
// The storage dnwp should be initialized as y0&y1 at the end.
I_CHANGESTAT_FN(i__intersect_net){
  ALLOC_AUX_STORAGE(1, StoreNetAndRefEL, storage);
  Network *dnwp = &(storage->nw);
  double *ref_el = storage->ref_el = INPUT_PARAM + 1;
  
  *dnwp = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, 0, 0, NULL);
  
  Edge nedges = *ref_el;
  for(Edge i=0; i<nedges; i++){
    Vertex tail=ref_el[1+i], head=ref_el[1+nedges+i];
    if(IS_OUTEDGE(tail, head)) {
      ToggleEdge(tail,head, dnwp);
    }
  }
}

U_CHANGESTAT_FN(u__intersect_net){
  GET_AUX_STORAGE(StoreNetAndRefEL, storage);
  Network *dnwp = &(storage->nw);
  double *ref_el = storage->ref_el;
  // only toggle if the edge is in y0. otherwise changing y1 won't matter.
  if(dEdgeListSearch(tail, head, ref_el))
    ToggleEdge(tail, head, dnwp);
}

F_CHANGESTAT_FN(f__intersect_net){
  GET_AUX_STORAGE(StoreNetAndRefEL, storage);
  Network *dnwp = &(storage->nw);
  
  NetworkDestroy(dnwp);
}

I_CHANGESTAT_FN(i__intersect_net_toggles_in_list){
  //Rprintf("allocating intersect_net_tog\n");
  ALLOC_AUX_STORAGE(1, StoreNetAndRefEL, storage);
  Network *dnwp = &(storage->nw);
  double *ref_el = storage->ref_el = INPUT_PARAM + 1;
  
  *dnwp = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, 0, 0, NULL);
  
  Edge nedges = *ref_el;
  for(Edge i=0; i<nedges; i++){
    Vertex tail=ref_el[1+i], head=ref_el[1+nedges+i];
    if(IS_OUTEDGE(tail, head)) {
      ToggleEdge(tail,head, dnwp);
    }
  }
}

U_CHANGESTAT_FN(u__intersect_net_toggles_in_list){
  GET_AUX_STORAGE(StoreNetAndRefEL, storage);
  Network *dnwp = &(storage->nw);
  
  ToggleEdge(tail, head, dnwp);
}

F_CHANGESTAT_FN(f__intersect_net_toggles_in_list){
  GET_AUX_STORAGE(StoreNetAndRefEL, storage);
  Network *dnwp = &(storage->nw);
  
  NetworkDestroy(dnwp);
}

// Initialize aux network to y1. Then loop through the edges of y0 (ref_el).
// If the edge does not exists in y1, then toggle it on in aux network.
// The storage dnwp should be initialized as y0|y1 at the end.
I_CHANGESTAT_FN(i__union_net){
  ALLOC_AUX_STORAGE(1, StoreNetAndRefEL, storage);
  Network *dnwp = &(storage->nw);
  double *ref_el = storage->ref_el = INPUT_PARAM + 1;
  
  NetworkCopy(dnwp, nwp);
  
  Edge nedges = *ref_el;
  for(Edge i=0; i<nedges; i++){
    Vertex tail=ref_el[1+i], head=ref_el[1+nedges+i];
    if(IS_OUTEDGE(tail, head)==0) {
      ToggleEdge(tail,head, dnwp);
    }
  }
}

U_CHANGESTAT_FN(u__union_net){
  GET_AUX_STORAGE(StoreNetAndRefEL, storage);
  Network *dnwp = &(storage->nw);
  double *ref_el = storage->ref_el;
  // If the edge is in y0, changing y1 won't matter.
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
