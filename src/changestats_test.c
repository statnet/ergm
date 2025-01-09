/*  File src/changestats_test.c in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_storage.h"
#include "ergm_changestats_auxnet.h"
#include "ergm_dyad_hashmap.h"

I_CHANGESTAT_FN(i_test_abs_edges_minus_5){
  ALLOC_STORAGE(1, Edge, edges);

  *edges = N_EDGES; // Pretend this takes a long time.
}

C_CHANGESTAT_FN(c_test_abs_edges_minus_5){
  GET_STORAGE(Edge, stored_edges_ptr);
  long int edges = *stored_edges_ptr;
  CHANGE_STAT[0] = -labs(edges-5);
  CHANGE_STAT[0] += labs(edges-5 + (edgestate?-1:1));
}

U_CHANGESTAT_FN(u_test_abs_edges_minus_5){
  GET_STORAGE(Edge, edges);
  *edges += edgestate ? - 1 : 1;
}

// Finalizer is not needed.


S_CHANGESTAT_FN(s_test_abs_edges_minus_5){
  GET_STORAGE(Edge, edges);
  // Storage uninitialized: compute from scratch.
  if(!edges){
    CHANGE_STAT[0] = labs((long int)N_EDGES-5);
  }else{ // Storage initialized: use it.
    CHANGE_STAT[0] = labs((long int)*edges-5);
  }
}

C_CHANGESTAT_FN(c_test_abs_edges_minus_5_no_s){c_test_abs_edges_minus_5(tail, head, mtp, nwp, edgestate);}
I_CHANGESTAT_FN(i_test_abs_edges_minus_5_no_s){i_test_abs_edges_minus_5(mtp, nwp);}
U_CHANGESTAT_FN(u_test_abs_edges_minus_5_no_s){u_test_abs_edges_minus_5(tail, head, mtp, nwp, edgestate);}

C_CHANGESTAT_FN(c_isociomatrix){
  GET_AUX_STORAGE(int *, sm);
  
    Dyad pos = tail-1 + (head-1)*N_NODES;
    CHANGE_STAT[pos] += sm[tail][head]? -1 : +1;
}

C_CHANGESTAT_FN(c_discord_isociomatrix){
  GET_AUX_STORAGE(int *, sm);
  
    Dyad pos = tail-1 + (head-1)*N_NODES;
    CHANGE_STAT[pos] += sm[tail][head]? -1 : +1;
}

C_CHANGESTAT_FN(c_disc_inter_union_net_Network){
  GET_AUX_STORAGE_NUM(StoreAuxnet, dauxnet, 0);
  GET_AUX_STORAGE_NUM(StoreAuxnet, iauxnet, 1);
  GET_AUX_STORAGE_NUM(StoreAuxnet, uauxnet, 2);

  int refedge = dEdgeListSearch(tail, head, INPUT_PARAM)!=0;
  
  CHANGE_STAT[0] = edgestate!=refedge ? -1 : +1;
  CHANGE_STAT[1] = refedge ? (edgestate ? -1 : +1) : 0;
  CHANGE_STAT[2] = !refedge ? (edgestate ? -1 : +1) : 0;

  //TODO: Implement test using MAP_TOGGLE.

  CHANGE_STAT[3] = (EDGECOUNT(dauxnet->onwp)+CHANGE_STAT[0])*(EDGECOUNT(dauxnet->onwp)+CHANGE_STAT[0]) - EDGECOUNT(dauxnet->onwp)*EDGECOUNT(dauxnet->onwp);
  CHANGE_STAT[4] = (EDGECOUNT(iauxnet->onwp)+CHANGE_STAT[1])*(EDGECOUNT(iauxnet->onwp)+CHANGE_STAT[1]) - EDGECOUNT(iauxnet->onwp)*EDGECOUNT(iauxnet->onwp);
  CHANGE_STAT[5] = (EDGECOUNT(uauxnet->onwp)+CHANGE_STAT[2])*(EDGECOUNT(uauxnet->onwp)+CHANGE_STAT[2]) - EDGECOUNT(uauxnet->onwp)*EDGECOUNT(uauxnet->onwp);
}

C_CHANGESTAT_FN(c_disc_inter_union_net_DyadSet){
  GET_AUX_STORAGE_NUM(StoreStrictDyadSetAndRefEL, dstorage, 0);
  GET_AUX_STORAGE_NUM(StoreStrictDyadSetAndRefEL, istorage, 1);
  GET_AUX_STORAGE_NUM(StoreStrictDyadSetAndRefEL, ustorage, 2);

  int refedge = dEdgeListSearch(tail, head, INPUT_PARAM)!=0;
  
  CHANGE_STAT[0] = edgestate!=refedge ? -1 : +1;
  CHANGE_STAT[1] = refedge ? (edgestate ? -1 : +1) : 0;
  CHANGE_STAT[2] = !refedge ? (edgestate ? -1 : +1) : 0;

  CHANGE_STAT[3] = (dstorage->nwp->size+CHANGE_STAT[0])*(dstorage->nwp->size+CHANGE_STAT[0]) - dstorage->nwp->size*dstorage->nwp->size;
  CHANGE_STAT[4] = (istorage->nwp->size+CHANGE_STAT[1])*(istorage->nwp->size+CHANGE_STAT[1]) - istorage->nwp->size*istorage->nwp->size;
  CHANGE_STAT[5] = (ustorage->nwp->size+CHANGE_STAT[2])*(ustorage->nwp->size+CHANGE_STAT[2]) - ustorage->nwp->size*ustorage->nwp->size;
}

/*****************
 changestat: c__edges_times
*****************/
C_CHANGESTAT_FN(c__edges_times) {
  CHANGE_STAT[0] += edgestate ? - *INPUT_PARAM : *INPUT_PARAM;
}

S_CHANGESTAT_FN(s__edges_tests) {
  CHANGE_STAT[0] = N_EDGES * *INPUT_PARAM;
}
