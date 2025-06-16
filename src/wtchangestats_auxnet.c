/*  File src/wtchangestats_auxnet.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#include "ergm_wtchangestat_auxnet.h"
#include "ergm_wtchangestats_auxnet.h"
/* #include "ergm_dyad_hashmap.h" */
/* #include "ergm_dyad_hashmap_utils.h" */

/* // sets aux network to y0 XOR y1 */
/* WtI_CHANGESTAT_FN(i__Wtdiscord_net_Network){ */
/*   I_WtAUXNET(NetworkCopy(nwp)); */
/*   int *ref_el = IINPUT_PARAM; */

/*   Edge nedges = *ref_el; */
/*   for(Edge i=0; i<nedges; i++){ */
/*     Vertex tail=ref_el[1+i], head=ref_el[1+nedges+i]; */
/*     ToggleEdge(tail,head, auxnet->onwp); */
/*   } */
/* } */

/* WtU_CHANGESTAT_FN(u__Wtdiscord_net_Network){ */
/*   GET_AUX_STORAGE(StoreAuxnet, auxnet); */
/*   ToggleEdge(tail, head, auxnet->onwp); */
/* } */

/* F_CHANGESTAT_FN(f__Wtdiscord_net_Network){ */
/*   GET_AUX_STORAGE(StoreAuxnet, auxnet); */
/*   NetworkDestroy(auxnet->onwp); */
/* } */

/* // Initialize empty aux network. Then loop through the edges of y0 (ref_el). */
/* // If the edge also exists in y1, then toggle it off in auxnet->onwp. */
/* // The storage auxnet->onwp should be initialized as y0&y1 at the end. */
/* WtI_CHANGESTAT_FN(i__Wtintersect_net_Network){ */
/*   I_WtAUXNET(NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE)); */
/*   int *ref_el = IINPUT_PARAM; */

/*   Edge nedges = *ref_el; */
/*   for(Edge i=0; i<nedges; i++){ */
/*     Vertex tail=ref_el[1+i], head=ref_el[1+nedges+i]; */
/*     if(IS_OUTEDGE(tail, head)) { */
/*       ToggleKnownEdge(tail,head, auxnet->onwp, FALSE); */
/*     } */
/*   } */
/* } */

/* WtU_CHANGESTAT_FN(u__Wtintersect_net_Network){ */
/*   GET_AUX_STORAGE(StoreAuxnet, auxnet); */
/*   int *ref_el = IINPUT_PARAM; */
/*   // only toggle if the edge is in y0. otherwise changing y1 won't matter. */
/*   if(iEdgeListSearch(tail, head, ref_el)) */
/*     ToggleEdge(tail, head, auxnet->onwp); */
/* } */

/* WtF_CHANGESTAT_FN(f__Wtintersect_net_Network){ */
/*   GET_AUX_STORAGE(StoreAuxnet, auxnet); */
/*   NetworkDestroy(auxnet->onwp); */
/* } */

/* WtI_CHANGESTAT_FN(i__Wtintersect_net_toggles_in_list_Network){ */
/*   //Rprintf("allocating intersect_net_tog\n"); */
/*   I_WtAUXNET(NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE)); */
/*   int *ref_el = IINPUT_PARAM; */

/*   Edge nedges = *ref_el; */
/*   for(Edge i=0; i<nedges; i++){ */
/*     Vertex tail=ref_el[1+i], head=ref_el[1+nedges+i]; */
/*     if(IS_OUTEDGE(tail, head)) { */
/*       ToggleEdge(tail,head, auxnet->onwp); */
/*     } */
/*   } */
/* } */

/* WtU_CHANGESTAT_FN(u__Wtintersect_net_toggles_in_list_Network){ */
/*   GET_AUX_STORAGE(StoreAuxnet, auxnet); */
/*   ToggleEdge(tail, head, auxnet->onwp); */
/* } */

/* WtF_CHANGESTAT_FN(f__Wtintersect_net_toggles_in_list_Network){ */
/*   GET_AUX_STORAGE(StoreAuxnet, auxnet); */
/*   NetworkDestroy(auxnet->onwp); */
/* } */

/* // Initialize aux network to y1. Then loop through the edges of y0 (ref_el). */
/* // If the edge does not exists in y1, then toggle it on in aux network. */
/* // The storage auxnet->onwp should be initialized as y0|y1 at the end. */
/* WtI_CHANGESTAT_FN(i__Wtunion_net_Network){ */
/*   I_WtAUXNET(NetworkCopy(nwp)); */
/*   int *ref_el = IINPUT_PARAM; */

/*   Edge nedges = *ref_el; */
/*   for(Edge i=0; i<nedges; i++){ */
/*     Vertex tail=ref_el[1+i], head=ref_el[1+nedges+i]; */
/*     if(IS_OUTEDGE(tail, head)==0) { */
/*       ToggleKnownEdge(tail,head, auxnet->onwp, FALSE); */
/*     } */
/*   } */
/* } */

/* WtU_CHANGESTAT_FN(u__Wtunion_net_Network){ */
/*   GET_AUX_STORAGE(StoreAuxnet, auxnet); */
/*   int *ref_el = IINPUT_PARAM; */
/*   // If the edge is in y0, changing y1 won't matter. */
/*   if(iEdgeListSearch(tail, head, ref_el)==0) */
/*     ToggleEdge(tail, head, auxnet->onwp); */
/* } */

/* WtF_CHANGESTAT_FN(f__Wtunion_net_Network){ */
/*   GET_AUX_STORAGE(StoreAuxnet, auxnet); */
/*   NetworkDestroy(auxnet->onwp); */
/* } */

/* /\* blockdiag_net: */
/*    maintains a network that mirrors the main network but excludes ties */
/*    not within specified blocks; also exports a numeric vector of block */
/*    memberships */
/* *\/ */

/* WtI_CHANGESTAT_FN(i__Wtblockdiag_net){ */
/*   I_WtAUXNET(NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE)); */
/*   int *b = IINPUT_PARAM-1; // tail and head are indexed from 1. */

/*   for(Vertex tail=1; tail <= N_TAILS; tail++){ */
/*     Vertex head; */
/*     Edge e; */
/*     STEP_THROUGH_OUTEDGES(tail, e, head) { */
/*       if(b[tail]==b[head]) */
/* 	ToggleKnownEdge(tail, head, auxnet->onwp, FALSE); */
/*     } */
/*   } */
/* } */

/* WtU_CHANGESTAT_FN(u__Wtblockdiag_net){ */
/*   GET_AUX_STORAGE(StoreAuxnet, auxnet); */
/*   int *b = IINPUT_PARAM-1; // tail and head are indexed from 1. */

/*   if(b[tail]==b[head]) */
/*     ToggleKnownEdge(tail, head, auxnet->onwp, edgestate); */
/* } */

/* WtF_CHANGESTAT_FN(f__Wtblockdiag_net){ */
/*   GET_AUX_STORAGE(StoreAuxnet, auxnet); */
/*   NetworkDestroy(auxnet->onwp); */
/* } */

/* _undir_net

   Maintain an undirected binary network symmetrized from LHS according to the specified rule:

   1: weak (max)
   2: strong (min)
   3: upper (tail < head)
   4: lower (tail > head)

*/

WtI_CHANGESTAT_FN(i__Wtundir_net){
  I_WtAUXNET(WtNetworkInitialize(NULL, NULL, NULL, 0, N_NODES, FALSE, BIPARTITE));

  unsigned int rule = IINPUT_PARAM[0];
  WtEXEC_THROUGH_NET_EDGES_PRE(tail, head, e, weight, {
      __Wtundir_net_totoggle;
      if(totoggle) WtSetEdge(tail, head, w, auxnet->onwp);
    });
}

WtU_CHANGESTAT_FN(u__Wtundir_net){
  GET_AUX_STORAGE(StoreWtAuxnet, auxnet);
  unsigned int rule = IINPUT_PARAM[0];

  __Wtundir_net_totoggle;
  if(totoggle) WtSetEdge(MIN(tail,head), MAX(tail,head),w, auxnet->onwp);
}

WtF_CHANGESTAT_FN(f__Wtundir_net){
  GET_AUX_STORAGE(StoreWtAuxnet, auxnet);
  WtNetworkDestroy(auxnet->onwp);
}

/* /\* _filter_formula_net  */

/*    Maintain a binary network that has an edge wherever the */
/*    contribution of the a given term (edges, nodematch, etc.) whose */
/*    dyadwise value is either 0 or 1 is 1. */

/* *\/ */

/* WtI_CHANGESTAT_FN(i__Wtfilter_formula_net){ */
/*   I_WtAUXNET(WtNetworkInitialize(NULL, NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE)); */
/*   GET_STORAGE(WtModel, m); */
/*   unsigned int op = IINPUT_PARAM[0]; */

/*   STORAGE = m = WtModelInitialize(getListElement(mtp->R, "submodel"), NULL, nwp, FALSE); */

/*   WtEXEC_THROUGH_NET_EDGES_PRE(t, h, w, e, { */
/*       Rboolean edgestate = TRUE; // We know the edge is present in nwp. */
/*       __Wtfilter_formula_net_totoggle(t, h, w, nwp, m); */
/*       if(totoggle) WtAddEdgeToTrees(t, h, w, auxnet->onwp); */
/*     }); */
/* } */

/* WtU_CHANGESTAT_FN(u__Wtfilter_formula_net){ */
/*   GET_AUX_STORAGE(StoreAuxnet, auxnet); */
/*   GET_STORAGE(Model, m); */
/*   unsigned int op = IINPUT_PARAM[0]; */

/*   __Wtfilter_formula_net_totoggle(tail, head, nwp, m); */
/*   if(totoggle){ */
/*     if(edgestate) WtDeleteEdgeFromTrees(tail,head,auxnet->onwp); */
/*     else WtAddEdgeToTrees(tail,head,weight,auxnet->onwp); */
/*   } */
/* } */

/* WtF_CHANGESTAT_FN(f__Wtfilter_formula_net){ */
/*   GET_AUX_STORAGE(StoreAuxnet, auxnet); */
/*   Model *m = STORAGE; */

/*   ModelDestroy(nwp, m); */
/*   STORAGE = NULL; */
/*   NetworkDestroy(auxnet->onwp); */
/*   // DestroyStats() will deallocate the rest. */
/* } */

/* _subgraph_*_net

   Maintain an induced subgraph with specified vertex subset (two
   subsets, if bipartite).

*/

WtI_CHANGESTAT_FN(i__Wtsubgraph_net){
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
    error("Error in i__Wtsubgraph_net(): unrecognised output network type.");
    break;
  }

  I_WtAUXNET(WtNetworkInitialize(NULL, NULL, NULL, 0, n, dir, bip));

  WtEXEC_THROUGH_NET_EDGES_PRE(tail, head, e, weight, {
      Vertex st = thmap[0][tail];
      Vertex sh = thmap[1][head];
      if(!DIRECTED & (st==0 || sh==0)){
        st = thmap[0][head];
        sh = thmap[1][tail];
      }
      if(st!=0 && sh!=0) WtAddEdgeToTrees(st, sh, weight, auxnet->onwp);
    });
}

WtU_CHANGESTAT_FN(u__Wtsubgraph_net){
  GET_AUX_STORAGE(StoreWtAuxnet, auxnet);
  GET_STORAGE(int*, thmap);
  Vertex st = thmap[0][tail];
  Vertex sh = thmap[1][head];
  if(!DIRECTED & (st==0 || sh==0)){
    st = thmap[0][head];
    sh = thmap[1][tail];
  }
  if(st!=0 && sh!=0) WtSetEdge(st, sh, weight, auxnet->onwp);
}

WtF_CHANGESTAT_FN(f__Wtsubgraph_net){
  GET_AUX_STORAGE(StoreWtAuxnet, auxnet);
  WtNetworkDestroy(auxnet->onwp);
  // DestroyStats() will deallocate the rest.
}

/* _transformed_net

   Maintain a valued network subject to the specified transformation

   1: square root

*/

WtI_CHANGESTAT_FN(i__Wttransformed_net){
  I_WtAUXNET(WtNetworkInitialize(NULL, NULL, NULL, 0, N_NODES, FALSE, BIPARTITE));

  unsigned int op = IINPUT_PARAM[0];
  WtEXEC_THROUGH_NET_EDGES_PRE(tail, head, e, weight, {
      __Wttransformed_net_totoggle;
      if(totoggle) WtSetEdge(tail, head, w, auxnet->onwp);
    });
}

WtU_CHANGESTAT_FN(u__Wttransformed_net){
  GET_AUX_STORAGE(StoreWtAuxnet, auxnet);
  unsigned int op = IINPUT_PARAM[0];

  __Wttransformed_net_totoggle;
  if(totoggle) WtSetEdge(tail, head, w, auxnet->onwp);
}

WtF_CHANGESTAT_FN(f__Wttransformed_net){
  GET_AUX_STORAGE(StoreWtAuxnet, auxnet);
  WtNetworkDestroy(auxnet->onwp);
}
