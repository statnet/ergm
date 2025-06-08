/*  File src/changestats_auxnet.c in package ergm, part of the Statnet suite of
 *  packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#include "ergm_changestat_auxnet.h"
#include "ergm_changestats_auxnet.h"
#include "ergm_dyad_hashmap.h"
#include "ergm_dyad_hashmap_utils.h"

// sets aux network to y0 XOR y1
I_CHANGESTAT_FN(i__discord_net_Network){
  I_AUXNET(NetworkCopy(nwp));
  int *ref_el = IINPUT_PARAM;

  Edge nedges = *ref_el;
  for(Edge i=0; i<nedges; i++){
    Vertex tail=ref_el[1+i], head=ref_el[1+nedges+i];
    ToggleEdge(tail,head, auxnet->onwp);
  }
}

U_CHANGESTAT_FN(u__discord_net_Network){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  ToggleEdge(tail, head, auxnet->onwp);
}

F_CHANGESTAT_FN(f__discord_net_Network){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  NetworkDestroy(auxnet->onwp);
}

// Initialize empty aux network. Then loop through the edges of y0 (ref_el).
// If the edge also exists in y1, then toggle it off in auxnet->onwp.
// The storage auxnet->onwp should be initialized as y0&y1 at the end.
I_CHANGESTAT_FN(i__intersect_net_Network){
  I_AUXNET(NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, 0, 0, NULL));
  int *ref_el = IINPUT_PARAM;

  Edge nedges = *ref_el;
  for(Edge i=0; i<nedges; i++){
    Vertex tail=ref_el[1+i], head=ref_el[1+nedges+i];
    if(IS_OUTEDGE(tail, head)) {
      ToggleKnownEdge(tail,head, auxnet->onwp, FALSE);
    }
  }
}

U_CHANGESTAT_FN(u__intersect_net_Network){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  int *ref_el = IINPUT_PARAM;
  // only toggle if the edge is in y0. otherwise changing y1 won't matter.
  if(iEdgeListSearch(tail, head, ref_el))
    ToggleEdge(tail, head, auxnet->onwp);
}

F_CHANGESTAT_FN(f__intersect_net_Network){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  NetworkDestroy(auxnet->onwp);
}

I_CHANGESTAT_FN(i__intersect_net_toggles_in_list_Network){
  //Rprintf("allocating intersect_net_tog\n");
  I_AUXNET(NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, 0, 0, NULL));
  int *ref_el = IINPUT_PARAM;

  Edge nedges = *ref_el;
  for(Edge i=0; i<nedges; i++){
    Vertex tail=ref_el[1+i], head=ref_el[1+nedges+i];
    if(IS_OUTEDGE(tail, head)) {
      ToggleEdge(tail,head, auxnet->onwp);
    }
  }
}

U_CHANGESTAT_FN(u__intersect_net_toggles_in_list_Network){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  ToggleEdge(tail, head, auxnet->onwp);
}

F_CHANGESTAT_FN(f__intersect_net_toggles_in_list_Network){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  NetworkDestroy(auxnet->onwp);
}

// Initialize aux network to y1. Then loop through the edges of y0 (ref_el).
// If the edge does not exists in y1, then toggle it on in aux network.
// The storage auxnet->onwp should be initialized as y0|y1 at the end.
I_CHANGESTAT_FN(i__union_net_Network){
  I_AUXNET(NetworkCopy(nwp));
  int *ref_el = IINPUT_PARAM;

  Edge nedges = *ref_el;
  for(Edge i=0; i<nedges; i++){
    Vertex tail=ref_el[1+i], head=ref_el[1+nedges+i];
    if(IS_OUTEDGE(tail, head)==0) {
      ToggleKnownEdge(tail,head, auxnet->onwp, FALSE);
    }
  }
}

U_CHANGESTAT_FN(u__union_net_Network){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  int *ref_el = IINPUT_PARAM;
  // If the edge is in y0, changing y1 won't matter.
  if(iEdgeListSearch(tail, head, ref_el)==0)
    ToggleEdge(tail, head, auxnet->onwp);
}

F_CHANGESTAT_FN(f__union_net_Network){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  NetworkDestroy(auxnet->onwp);
}

/* blockdiag_net:
   maintains a network that mirrors the main network but excludes ties
   not within specified blocks; also exports a numeric vector of block
   memberships
*/

I_CHANGESTAT_FN(i__blockdiag_net){
  I_AUXNET(NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, 0, 0, NULL));
  int *b = IINPUT_PARAM-1; // tail and head are indexed from 1.

  for(Vertex tail=1; tail <= N_TAILS; tail++){
    Vertex head;
    Edge e;
    STEP_THROUGH_OUTEDGES(tail, e, head) {
      if(b[tail]==b[head])
	ToggleKnownEdge(tail, head, auxnet->onwp, FALSE);
    }
  }
}

U_CHANGESTAT_FN(u__blockdiag_net){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  int *b = IINPUT_PARAM-1; // tail and head are indexed from 1.

  if(b[tail]==b[head])
    ToggleKnownEdge(tail, head, auxnet->onwp, edgestate);
}

F_CHANGESTAT_FN(f__blockdiag_net){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  NetworkDestroy(auxnet->onwp);
}

/* _undir_net

   Maintain an undirected binary network symmetrized from LHS according to the specified rule:

   1: weak
   2: strong
   3: upper (tail < head)
   4: lower (tail > head)

*/

I_CHANGESTAT_FN(i__undir_net){
  I_AUXNET(NetworkInitialize(NULL, NULL, 0, N_NODES, FALSE, BIPARTITE, FALSE, 0, NULL));

  unsigned int rule = IINPUT_PARAM[0];
  EXEC_THROUGH_NET_EDGES_PRE(tail, head, e, {
      __undir_net_totoggle;
      if(totoggle && !IS_UNDIRECTED_EDGE(tail,head,auxnet->onwp)) ToggleKnownEdge(tail, head, auxnet->onwp, FALSE);
    });
}

U_CHANGESTAT_FN(u__undir_net){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  unsigned int rule = IINPUT_PARAM[0];

  __undir_net_totoggle;
  if(totoggle) ToggleEdge(MIN(tail,head),MAX(tail,head),auxnet->onwp);
}

F_CHANGESTAT_FN(f__undir_net){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  NetworkDestroy(auxnet->onwp);
}

/* _filter_formula_net

   Maintain a binary network that has an edge wherever the
   contribution of the a given term (edges, nodematch, etc.) whose
   dyadwise value is either 0 or 1 is 1.

*/

I_CHANGESTAT_FN(i__filter_formula_net){
  I_AUXNET(NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, FALSE, 0, NULL));
  GET_STORAGE(Model, m);
  unsigned int op = IINPUT_PARAM[0];

  STORAGE = m = ModelInitialize(getListElement(mtp->R, "submodel"), NULL, nwp, FALSE);

  EXEC_THROUGH_NET_EDGES_PRE(t, h, e, {
      Rboolean edgestate = TRUE; // We know the edge is present in nwp.
      __filter_formula_net_totoggle(t, h, nwp, m);
      if(totoggle) AddEdgeToTrees(t, h, auxnet->onwp);
    });
}

U_CHANGESTAT_FN(u__filter_formula_net){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  GET_STORAGE(Model, m);
  unsigned int op = IINPUT_PARAM[0];

  __filter_formula_net_totoggle(tail, head, nwp, m);
  if(totoggle){
    if(edgestate) DeleteEdgeFromTrees(tail,head,auxnet->onwp);
    else AddEdgeToTrees(tail,head,auxnet->onwp);
  }
}

F_CHANGESTAT_FN(f__filter_formula_net){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  Model *m = STORAGE;

  ModelDestroy(nwp, m);
  STORAGE = NULL;
  NetworkDestroy(auxnet->onwp);
  // DestroyStats() will deallocate the rest.
}

/* _subgraph_*_net

   Maintain an induced subgraph with specified vertex subset (two
   subsets, if bipartite).

*/

I_CHANGESTAT_FN(i__subgraph_net){
  ALLOC_STORAGE(2, int*, thmap);
  int *inputs = IINPUT_PARAM;
  unsigned int type = *(inputs++);
  /* These will be overwritten in the following switch statement. */
  Vertex n=0, bip=0;
  Rboolean dir=FALSE;

  switch(type){
  case 1:
    bip = 0;
    n = *(inputs++);
    dir = TRUE;
    thmap[0] = thmap[1] = inputs - 1;
    break;
  case 2:
    bip = 0;
    n = *(inputs++);
    dir = FALSE;
    thmap[0] = thmap[1] = inputs - 1;
    break;
  case 3:
    bip = *(inputs++);
    n = bip + *(inputs++);
    thmap[0] = inputs - 1;
    thmap[1] = inputs + nwp->nnodes - 1;
    break;
  default:
    error("Error in i__subgraph_net(): unrecognised output network type.");
    break;
  }

  I_AUXNET(NetworkInitialize(NULL, NULL, 0, n, dir, bip, FALSE, 0, NULL));

  EXEC_THROUGH_NET_EDGES_PRE(tail, head, e, {
      Vertex st = thmap[0][tail];
      Vertex sh = thmap[1][head];
      if(!DIRECTED & (st==0 || sh==0)){
        st = thmap[0][head];
        sh = thmap[1][tail];
      }
      if(st!=0 && sh!=0) AddEdgeToTrees(st, sh, auxnet->onwp);
    });
}

U_CHANGESTAT_FN(u__subgraph_net){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  GET_STORAGE(int*, thmap);
  Vertex st = thmap[0][tail];
  Vertex sh = thmap[1][head];
  if(!DIRECTED & (st==0 || sh==0)){
    st = thmap[0][head];
    sh = thmap[1][tail];
  }
  if(st!=0 && sh!=0) ToggleKnownEdge(st, sh, auxnet->onwp, edgestate);
}

F_CHANGESTAT_FN(f__subgraph_net){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  NetworkDestroy(auxnet->onwp);
  // DestroyStats() will deallocate the rest.
}
