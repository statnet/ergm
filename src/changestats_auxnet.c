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
  GET_AUX_STORAGE_NUM(StoreNetAndRefEL, storage, 1);

  nwp = storage->nwp; // So that we can use the macros.
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

// sets aux network to y0 XOR y1
I_CHANGESTAT_FN(i__discord_net_Network){
  ALLOC_AUX_STORAGE(1, StoreNetAndRefEL, storage);
  Network *dnwp = storage->nwp = NetworkCopy(nwp);
  double *ref_el = storage->ref_el = INPUT_PARAM + 1;
  
  Edge nedges = *ref_el;
  for(Edge i=0; i<nedges; i++){
    Vertex tail=ref_el[1+i], head=ref_el[1+nedges+i];
    ToggleEdge(tail,head, dnwp);
  }
}

U_CHANGESTAT_FN(u__discord_net_Network){
  GET_AUX_STORAGE(StoreNetAndRefEL, storage);
  Network *dnwp = storage->nwp;

  ToggleEdge(tail, head, dnwp);
}

F_CHANGESTAT_FN(f__discord_net_Network){
  GET_AUX_STORAGE(StoreNetAndRefEL, storage);
  Network *dnwp = storage->nwp;

  NetworkDestroy(dnwp);
}

// Initialize empty aux network. Then loop through the edges of y0 (ref_el).
// If the edge also exists in y1, then toggle it off in dnwp.
// The storage dnwp should be initialized as y0&y1 at the end.
I_CHANGESTAT_FN(i__intersect_net_Network){
  ALLOC_AUX_STORAGE(1, StoreNetAndRefEL, storage);
  Network *dnwp = storage->nwp = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, 0, 0, NULL);
  double *ref_el = storage->ref_el = INPUT_PARAM + 1;
  
  Edge nedges = *ref_el;
  for(Edge i=0; i<nedges; i++){
    Vertex tail=ref_el[1+i], head=ref_el[1+nedges+i];
    if(IS_OUTEDGE(tail, head)) {
      ToggleEdge(tail,head, dnwp);
    }
  }
}

U_CHANGESTAT_FN(u__intersect_net_Network){
  GET_AUX_STORAGE(StoreNetAndRefEL, storage);
  Network *dnwp = storage->nwp;
  double *ref_el = storage->ref_el;
  // only toggle if the edge is in y0. otherwise changing y1 won't matter.
  if(dEdgeListSearch(tail, head, ref_el))
    ToggleEdge(tail, head, dnwp);
}

F_CHANGESTAT_FN(f__intersect_net_Network){
  GET_AUX_STORAGE(StoreNetAndRefEL, storage);
  Network *dnwp = storage->nwp;

  NetworkDestroy(dnwp);
}

I_CHANGESTAT_FN(i__intersect_net_toggles_in_list_Network){
  //Rprintf("allocating intersect_net_tog\n");
  ALLOC_AUX_STORAGE(1, StoreNetAndRefEL, storage);
  Network *dnwp = storage->nwp = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, 0, 0, NULL);
  double *ref_el = storage->ref_el = INPUT_PARAM + 1;
  
  Edge nedges = *ref_el;
  for(Edge i=0; i<nedges; i++){
    Vertex tail=ref_el[1+i], head=ref_el[1+nedges+i];
    if(IS_OUTEDGE(tail, head)) {
      ToggleEdge(tail,head, dnwp);
    }
  }
}

U_CHANGESTAT_FN(u__intersect_net_toggles_in_list_Network){
  GET_AUX_STORAGE(StoreNetAndRefEL, storage);
  Network *dnwp = storage->nwp;

  ToggleEdge(tail, head, dnwp);
}

F_CHANGESTAT_FN(f__intersect_net_toggles_in_list_Network){
  GET_AUX_STORAGE(StoreNetAndRefEL, storage);
  Network *dnwp = storage->nwp;

  NetworkDestroy(dnwp);
}

// Initialize aux network to y1. Then loop through the edges of y0 (ref_el).
// If the edge does not exists in y1, then toggle it on in aux network.
// The storage dnwp should be initialized as y0|y1 at the end.
I_CHANGESTAT_FN(i__union_net_Network){
  ALLOC_AUX_STORAGE(1, StoreNetAndRefEL, storage);
  Network *dnwp = storage->nwp = NetworkCopy(nwp);
  double *ref_el = storage->ref_el = INPUT_PARAM + 1;
  
  Edge nedges = *ref_el;
  for(Edge i=0; i<nedges; i++){
    Vertex tail=ref_el[1+i], head=ref_el[1+nedges+i];
    if(IS_OUTEDGE(tail, head)==0) {
      ToggleEdge(tail,head, dnwp);
    }
  }
}

U_CHANGESTAT_FN(u__union_net_Network){
  GET_AUX_STORAGE(StoreNetAndRefEL, storage);
  Network *dnwp = storage->nwp;
  double *ref_el = storage->ref_el;
  // If the edge is in y0, changing y1 won't matter.
  if(dEdgeListSearch(tail, head, ref_el)==0)
    ToggleEdge(tail, head, dnwp);
}

F_CHANGESTAT_FN(f__union_net_Network){
  GET_AUX_STORAGE(StoreNetAndRefEL, storage);
  Network *dnwp = storage->nwp;

  NetworkDestroy(dnwp);
}

I_CHANGESTAT_FN(i__discord_net_DyadSet){
  ALLOC_AUX_STORAGE(1, StoreDyadSetAndRefEL, storage);
  StoreDyadSet *dnwp = storage->nwp = NetworkToDyadSet(nwp);
  double *ref_el = storage->ref_el = INPUT_PARAM + 1;
  
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
  double *ref_el = storage->ref_el = INPUT_PARAM + 1;
  
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
  double *ref_el = storage->ref_el;
  // only toggle if the edge is in y0. otherwise changing y1 won't matter.
  if(dEdgeListSearch(tail, head, ref_el))
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
  double *ref_el = storage->ref_el = INPUT_PARAM + 1;
  
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
  double *ref_el = storage->ref_el = INPUT_PARAM + 1;
  
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
  double *ref_el = storage->ref_el;
  // If the edge is in y0, changing y1 won't matter.
  if(dEdgeListSearch(tail, head, ref_el)==0)
    DyadSetToggle(tail,head, dnwp);
}

F_CHANGESTAT_FN(f__union_net_DyadSet){
  GET_AUX_STORAGE(StoreDyadSetAndRefEL, storage);
  StoreDyadSet *dnwp = storage->nwp;

  kh_destroy(DyadSet, dnwp);
}

/* blockdiag_net:
   maintains a network that mirrors the main network but excludes ties
   not within specified blocks; also exports a numeric vector of block
   memberships
 */

I_CHANGESTAT_FN(i__blockdiag_net){
  ALLOC_AUX_STORAGE(1, StoreNetAndBID, blkinfo);
  blkinfo->nwp = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, 0, 0, NULL);
  Network *bnwp = blkinfo->nwp;
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
  Network *bnwp = blkinfo->nwp;
  double *b = blkinfo->b;

  if(b[tail]==b[head])
    ToggleEdge(tail, head, bnwp);
}

F_CHANGESTAT_FN(f__blockdiag_net){
  GET_AUX_STORAGE(StoreNetAndBID, blkinfo);
  Network *bnwp = blkinfo->nwp;
  NetworkDestroy(bnwp);
}

/* _undir_net 

   Maintain an undirected binary network symmetrized from LHS according to the specified rule:

   1: weak
   2: strong
   3: upper (tail < head)
   4: lower (tail > head)

*/
I_CHANGESTAT_FN(i__undir_net){
  double *inputs = INPUT_PARAM;
  Network *unwp = AUX_STORAGE = NetworkInitialize(NULL, NULL, 0, N_NODES, FALSE, BIPARTITE, FALSE, 0, NULL); inputs++;
  unsigned int rule = *inputs;
   
  EXEC_THROUGH_NET_EDGES_PRE(t, h, e, {
      unsigned int toadd;
      switch(rule){
      case 1: // weak
	toadd = TRUE;
	break;
      case 2: // strong
	toadd = IS_OUTEDGE(h,t);
	break;
      case 3: // upper
	toadd = t<=h;
	break;
      case 4: // lower
	toadd = t>=h;
	break;
      default: // never reached, but avoids a warning
	toadd = FALSE;
      }
      if(toadd) AddEdgeToTrees(MIN(t,h), MAX(t,h), unwp);
    });
}

U_CHANGESTAT_FN(u__undir_net){
  double *inputs = INPUT_PARAM;
  GET_AUX_STORAGE(Network, unwp); inputs++;
  unsigned int rule = *inputs;

  unsigned int totoggle;
  switch(rule){
  case 1: // weak
    totoggle = !IS_OUTEDGE(head,tail);
    break;
  case 2: // strong
    totoggle = IS_OUTEDGE(head,tail);
    break;
  case 3: // upper
    totoggle = tail<=head;
    break;
  case 4: // lower
    totoggle = tail>=head;
    break;
  default: // never reached, but avoids a warning
    totoggle = FALSE;
  }
  if(totoggle) ToggleEdge(MIN(tail,head),MAX(tail,head),unwp);
}

F_CHANGESTAT_FN(f__undir_net){
  GET_AUX_STORAGE(Network, unwp);
  NetworkDestroy(unwp);
  AUX_STORAGE=NULL;
  // DestroyStats() will deallocate the rest.
}

/* _subgraph_net

   Maintain an induced subgraph with specified vertex subset (two
   subsets, if bipartite).

*/
I_CHANGESTAT_FN(i__subgraph_net){
  double *inputs = INPUT_PARAM;
  ALLOC_AUX_STORAGE(1, StoreSubgraph, storage); inputs++;
  unsigned int type = *(inputs++);
  switch(type){
  case 1: // directed
    storage->nwp = NetworkInitialize(NULL, NULL, 0, *inputs, TRUE, FALSE, FALSE, 0, NULL);
    inputs++;
    storage->tmap = storage->hmap = inputs - 1;
    break;
  case 2: // undirected
    storage->nwp = NetworkInitialize(NULL, NULL, 0, *inputs, FALSE, FALSE, FALSE, 0, NULL);
    inputs++;
    storage->tmap = storage->hmap = inputs - 1;
    break;
  case 3: // bipartite undirected
    storage->nwp = NetworkInitialize(NULL, NULL, 0, inputs[0]+inputs[1], FALSE, inputs[0], FALSE, 0, NULL);
    inputs+=2;
    storage->tmap = inputs - 1;
    storage->hmap = inputs + nwp->nnodes - 1;
    break;
  default: // never reached, but avoids a warning
    break;
  }
  EXEC_THROUGH_NET_EDGES_PRE(t, h, e, {
      Vertex st = storage->tmap[t];
      Vertex sh = storage->hmap[h];
      if(!DIRECTED && (st==0 || sh==0)){
	st = storage->tmap[h];
	sh = storage->hmap[t];
      }
      if(st!=0 && sh!=0) AddEdgeToTrees(st, sh, storage->nwp);
    });
}

U_CHANGESTAT_FN(u__subgraph_net){
  GET_AUX_STORAGE(StoreSubgraph, storage);
  Vertex st = storage->tmap[tail];
  Vertex sh = storage->hmap[head];
  if(!DIRECTED && (st==0 || sh==0)){
    st = storage->tmap[head];
    sh = storage->hmap[tail];
  }
  if(st!=0 && sh!=0) ToggleKnownEdge(st, sh, storage->nwp, edgeflag);
}

F_CHANGESTAT_FN(f__subgraph_net){
  GET_AUX_STORAGE(StoreSubgraph, storage);
  NetworkDestroy(storage->nwp);
  // DestroyStats() will deallocate the rest.
}
