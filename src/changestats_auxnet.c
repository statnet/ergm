/*  File src/changestats_test.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2017 Statnet Commons
 */
#include "ergm_changestat_auxnet.h"
#include "ergm_changestat_operator.h"

MAP_TOGGLE_FN(map_toggle__discord_net_Network){
  MAP_TOGGLE_MAXTOGGLES(1);
  MAP_TOGGLE_PROPAGATE;
}

// sets aux network to y0 XOR y1
I_CHANGESTAT_FN(i__discord_net_Network){
  I_AUXNET(NetworkCopy(nwp), map_toggle__discord_net_Network);
  double *ref_el = INPUT_PARAM + 1;
  
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


MAP_TOGGLE_FN(map_toggle__intersect_net_Network){
  MAP_TOGGLE_MAXTOGGLES(1);

  ModelTerm *mtp = auxnet->mtp;
  double *ref_el = INPUT_PARAM + 1;
  MAP_TOGGLE_PROPAGATE_IF(dEdgeListSearch(tail, head, ref_el));
}

// Initialize empty aux network. Then loop through the edges of y0 (ref_el).
// If the edge also exists in y1, then toggle it off in auxnet->onwp.
// The storage auxnet->onwp should be initialized as y0&y1 at the end.
I_CHANGESTAT_FN(i__intersect_net_Network){
  I_AUXNET(NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, 0, 0, NULL), map_toggle__intersect_net_Network);
  double *ref_el = INPUT_PARAM + 1;
  
  Edge nedges = *ref_el;
  for(Edge i=0; i<nedges; i++){
    Vertex tail=ref_el[1+i], head=ref_el[1+nedges+i];
    if(IS_OUTEDGE(tail, head)) {
      ToggleEdge(tail,head, auxnet->onwp);
    }
  }
}

U_CHANGESTAT_FN(u__intersect_net_Network){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  double *ref_el = INPUT_PARAM + 1;
  // only toggle if the edge is in y0. otherwise changing y1 won't matter.
  if(dEdgeListSearch(tail, head, ref_el))
    ToggleEdge(tail, head, auxnet->onwp);
}

F_CHANGESTAT_FN(f__intersect_net_Network){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  NetworkDestroy(auxnet->onwp);
}

I_CHANGESTAT_FN(i__intersect_net_toggles_in_list_Network){
  //Rprintf("allocating intersect_net_tog\n");
  I_AUXNET(NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, 0, 0, NULL), map_toggle__discord_net_Network);
  double *ref_el = INPUT_PARAM + 1;
  
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

MAP_TOGGLE_FN(map_toggle__union_net_Network){
  MAP_TOGGLE_MAXTOGGLES(1);

  ModelTerm *mtp = auxnet->mtp;
  double *ref_el = INPUT_PARAM + 1;
  MAP_TOGGLE_PROPAGATE_IF(!dEdgeListSearch(tail, head, ref_el));
}

// Initialize aux network to y1. Then loop through the edges of y0 (ref_el).
// If the edge does not exists in y1, then toggle it on in aux network.
// The storage auxnet->onwp should be initialized as y0|y1 at the end.
I_CHANGESTAT_FN(i__union_net_Network){
  I_AUXNET(NetworkCopy(nwp), map_toggle__union_net_Network); 
  double *ref_el = INPUT_PARAM + 1;
  
  Edge nedges = *ref_el;
  for(Edge i=0; i<nedges; i++){
    Vertex tail=ref_el[1+i], head=ref_el[1+nedges+i];
    if(IS_OUTEDGE(tail, head)==0) {
      ToggleEdge(tail,head, auxnet->onwp);
    }
  }
}

U_CHANGESTAT_FN(u__union_net_Network){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  double *ref_el = INPUT_PARAM + 1;
  // If the edge is in y0, changing y1 won't matter.
  if(dEdgeListSearch(tail, head, ref_el)==0)
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

MAP_TOGGLE_FN(map_toggle__blockdiag_net){
  MAP_TOGGLE_MAXTOGGLES(1);

  ModelTerm *mtp = auxnet->mtp;
  double *b = INPUT_PARAM;
  MAP_TOGGLE_PROPAGATE_IF(b[tail]==b[head]);
}

I_CHANGESTAT_FN(i__blockdiag_net){
  I_AUXNET(NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, 0, 0, NULL), map_toggle__blockdiag_net);
  double *b = INPUT_PARAM;

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
  double *b = INPUT_PARAM;

  if(b[tail]==b[head])
    ToggleKnownEdge(tail, head, auxnet->onwp, edgeflag);
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

#define _undir_net_totoggle				\
  Rboolean totoggle;					\
  switch(rule){						\
  case 1: /* weak */					\
    totoggle = !IS_OUTEDGE(head,tail);			\
    break;						\
  case 2: /* strong */					\
    totoggle = IS_OUTEDGE(head,tail);			\
    break;						\
  case 3: /* upper */					\
    totoggle = tail<=head;				\
    break;						\
  case 4: /* lower */					\
    totoggle = tail>=head;				\
    break;						\
  default: /* never reached, but avoids a warning */	\
    totoggle = FALSE;					\
  }

MAP_TOGGLE_FN(map_toggle__undir_net){
  MAP_TOGGLE_MAXTOGGLES(1);

  ModelTerm *mtp = auxnet->mtp;
  unsigned int rule = INPUT_PARAM[1];
  Network *nwp = auxnet->inwp;
  _undir_net_totoggle;  
  if(totoggle){
    *tails = MIN(tail,head);
    *heads = MAX(tail,head);
    return 1;
  }else{
    return 0;
  }
}

I_CHANGESTAT_FN(i__undir_net){
  I_AUXNET(NetworkInitialize(NULL, NULL, 0, N_NODES, FALSE, BIPARTITE, FALSE, 0, NULL), map_toggle__undir_net);

  unsigned int rule = INPUT_PARAM[1];
  EXEC_THROUGH_NET_EDGES_PRE(tail, head, e, {
      _undir_net_totoggle;
      if(totoggle) AddEdgeToTrees(MIN(tail,head), MAX(tail,head), auxnet->onwp);
    });
}

U_CHANGESTAT_FN(u__undir_net){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  unsigned int rule = INPUT_PARAM[1];

  _undir_net_totoggle;
  if(totoggle) ToggleEdge(MIN(tail,head),MAX(tail,head),auxnet->onwp);
}

F_CHANGESTAT_FN(f__undir_net){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  NetworkDestroy(auxnet->onwp);
}

#undef _undir_net_totoggle


/* _filter_formula_net 

   Maintain a binary network that has an edge wherever the
   contribution of the a given term (edges, nodematch, etc.) whose
   dyadwise value is either 0 or 1 is 1.

*/


MAP_TOGGLE_FN(map_toggle__filter_formula_net){
  MAP_TOGGLE_MAXTOGGLES(1);

  ModelTerm *mtp = auxnet->mtp;
  GET_STORAGE(Model, m);
  ChangeStats(1, &tail, &head, auxnet->inwp, m);
  MAP_TOGGLE_PROPAGATE_IF(*(m->workspace)!=0);
}

I_CHANGESTAT_FN(i__filter_formula_net){
  I_AUXNET(NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, FALSE, 0, NULL), map_toggle__filter_formula_net);
  GET_STORAGE(Model, m);

  double *inputs = INPUT_PARAM+1;
  STORAGE = m = unpack_Model_as_double(&inputs);

  InitStats(nwp, m);

  EXEC_THROUGH_NET_EDGES_PRE(t, h, e, {
      ChangeStats(1, &t, &h, nwp, m);
      // I.e., if toggling the dyad changes the statistic, add
      // edge to the filter network.
      if(*(m->workspace)!=0) 
	AddEdgeToTrees(t, h, auxnet->onwp);
    });
}

U_CHANGESTAT_FN(u__filter_formula_net){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  Model *m = STORAGE;

  ChangeStats(1, &tail, &head, nwp, m);
  if(*(m->workspace)!=0){
    if(edgeflag) DeleteEdgeFromTrees(tail,head,auxnet->onwp);
    else AddEdgeToTrees(tail,head,auxnet->onwp);
  }

  UPDATE_STORAGE(tail, head, nwp, m, NULL, edgeflag);
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

/* directed */

MAP_TOGGLE_FN(map_toggle__subgraph_dir_net){
  MAP_TOGGLE_MAXTOGGLES(1);
  ModelTerm *mtp = auxnet->mtp;
  GET_STORAGE(double*, thmap);

  Vertex st = thmap[0][tail];
  Vertex sh = thmap[1][head];
  if(st!=0 && sh!=0){
    *tails = st;
    *heads = sh;
    return 1;
  }else{
    return 0;
  }
}

I_CHANGESTAT_FN(i__subgraph_dir_net){
  ALLOC_STORAGE(2, double*, thmap);
  double *inputs = INPUT_PARAM+1;
  I_AUXNET(NetworkInitialize(NULL, NULL, 0, *inputs, TRUE, FALSE, FALSE, 0, NULL), map_toggle__subgraph_dir_net);
  inputs++;
  thmap[0] = thmap[1] = inputs - 1;

  EXEC_THROUGH_NET_EDGES_PRE(tail, head, e, {
      Vertex st = thmap[0][tail];
      Vertex sh = thmap[1][head];
      if(st!=0 && sh!=0) AddEdgeToTrees(st, sh, auxnet->onwp);
    });
}

U_CHANGESTAT_FN(u__subgraph_dir_net){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  GET_STORAGE(double*, thmap);
  Vertex st = thmap[0][tail];
  Vertex sh = thmap[1][head];
  if(st!=0 && sh!=0) ToggleKnownEdge(st, sh, auxnet->onwp, edgeflag);
}

F_CHANGESTAT_FN(f__subgraph_dir_net){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  NetworkDestroy(auxnet->onwp);
  // DestroyStats() will deallocate the rest.
}

/* undirected */

MAP_TOGGLE_FN(map_toggle__subgraph_undir_net){
  MAP_TOGGLE_MAXTOGGLES(1);
  ModelTerm *mtp = auxnet->mtp;
  GET_STORAGE(double*, thmap);

  Vertex st = thmap[0][tail];
  Vertex sh = thmap[1][head];
  if(st==0 || sh==0){
    st = thmap[0][head];
    sh = thmap[1][tail];
  }
  if(st!=0 && sh!=0){
    *tails = st;
    *heads = sh;
    return 1;
  }else{
    return 0;
  }
}

I_CHANGESTAT_FN(i__subgraph_undir_net){
  ALLOC_STORAGE(2, double*, thmap);
  double *inputs = INPUT_PARAM+1;
  I_AUXNET(NetworkInitialize(NULL, NULL, 0, *inputs, FALSE, FALSE, FALSE, 0, NULL), map_toggle__subgraph_undir_net);
  inputs++;
  thmap[0] = thmap[1] = inputs - 1;

  EXEC_THROUGH_NET_EDGES_PRE(tail, head, e, {
      Vertex st = thmap[0][tail];
      Vertex sh = thmap[1][head];
      if(st==0 || sh==0){
        st = thmap[0][head];
        sh = thmap[1][tail];
      }
      if(st!=0 && sh!=0) AddEdgeToTrees(st, sh, auxnet->onwp);
    });
}

U_CHANGESTAT_FN(u__subgraph_undir_net){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  GET_STORAGE(double*, thmap);
  Vertex st = thmap[0][tail];
  Vertex sh = thmap[1][head];
  if(st==0 || sh==0){
    st = thmap[0][head];
    sh = thmap[1][tail];
  }
  if(st!=0 && sh!=0) ToggleKnownEdge(st, sh, auxnet->onwp, edgeflag);
}

F_CHANGESTAT_FN(f__subgraph_undir_net){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  NetworkDestroy(auxnet->onwp);
  // DestroyStats() will deallocate the rest.
}

/* undirected bipartite */

MAP_TOGGLE_FN(map_toggle__subgraph_bip_net){
  MAP_TOGGLE_MAXTOGGLES(1);
  Network *nwp = auxnet->inwp;
  ModelTerm *mtp = auxnet->mtp;
  GET_STORAGE(double*, thmap);

  Vertex st = thmap[0][tail];
  Vertex sh = thmap[1][head];
  if(!DIRECTED && (st==0 || sh==0)){
    st = thmap[0][head];
    sh = thmap[1][tail];
  }
  if(st!=0 && sh!=0){
    *tails = st;
    *heads = sh;
    return 1;
  }else{
    return 0;
  }
}

I_CHANGESTAT_FN(i__subgraph_bip_net){
  ALLOC_STORAGE(2, double*, thmap);
  double *inputs = INPUT_PARAM+1;
  I_AUXNET(NetworkInitialize(NULL, NULL, 0, inputs[0]+inputs[1], FALSE, inputs[0], FALSE, 0, NULL), map_toggle__subgraph_bip_net);
  inputs+=2;
  thmap[0] = inputs - 1;
  thmap[1] = inputs + nwp->nnodes - 1;

  EXEC_THROUGH_NET_EDGES_PRE(tail, head, e, {
      Vertex st = thmap[0][tail];
      Vertex sh = thmap[1][head];
      if(st==0 || sh==0){
        st = thmap[0][head];
        sh = thmap[1][tail];
      }
      if(st!=0 && sh!=0) AddEdgeToTrees(st, sh, auxnet->onwp);
    });
}

U_CHANGESTAT_FN(u__subgraph_bip_net){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  GET_STORAGE(double*, thmap);
  Vertex st = thmap[0][tail];
  Vertex sh = thmap[1][head];
  if(st==0 || sh==0){
    st = thmap[0][head];
    sh = thmap[1][tail];
  }
  if(st!=0 && sh!=0) ToggleKnownEdge(st, sh, auxnet->onwp, edgeflag);
}

F_CHANGESTAT_FN(f__subgraph_bip_net){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  NetworkDestroy(auxnet->onwp);
  // DestroyStats() will deallocate the rest.
}
