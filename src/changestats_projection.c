/*  File src/changestats_projection.c in package ergm, part of the Statnet
 *  suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */
#define STRICT_Wt_HEADERS
#include "ergm_wtchangestat_operator.h"
#include "ergm_changestat_operator.h"
#include "ergm_storage.h"

typedef struct {
  WtModel *m;
  Vertex *t, *h;
  double *w;
} StoreWtModelAndWtChanges;


I_CHANGESTAT_FN(i_on_proj_net){
  GET_AUX_STORAGE(WtNetwork, pnwp);
  ALLOC_STORAGE(1, StoreWtModelAndWtChanges, storage);
  // No need to allocate it: we are only storing a pointer to a model.
  storage->m = WtModelInitialize(getListElement(mtp->R, "submodel"), mtp->ext_state, pnwp, FALSE);

  /* WtSELECT_C_OR_D_BASED_ON_SUBMODEL(storage->m); */
  WtDELETE_IF_UNUSED_IN_SUBMODEL(x_func, storage->m);
  WtDELETE_IF_UNUSED_IN_SUBMODEL(z_func, storage->m);

  unsigned int mode = IINPUT_PARAM[0],
    maxchanges = mode == 1 ? BIPARTITE-1 : N_NODES-BIPARTITE-1;

  storage->t = R_Calloc(maxchanges, Vertex);
  storage->h = R_Calloc(maxchanges, Vertex);
  storage->w = R_Calloc(maxchanges, double);
}

/* TODO: A d_function? */

C_CHANGESTAT_FN(c_on_proj_net){
  GET_AUX_STORAGE(WtNetwork, pnwp);
  GET_STORAGE(StoreWtModelAndWtChanges, storage);

  int echange = edgestate ? -1 : 1;

  unsigned int nt = 0, mode = IINPUT_PARAM[0];

  switch(mode){
  case 1:
    EXEC_THROUGH_FINEDGES(head, e, tail2, {
        if(tail!=tail2){
          storage->t[nt] = MIN(tail, tail2);
          storage->h[nt] = MAX(tail, tail2);
          storage->w[nt] = WtGETWT(tail, tail2, pnwp) + echange;
          nt++;
        }
      }); break;
  case 2:
    EXEC_THROUGH_FOUTEDGES(tail, e, head2, {
        if(head!=head2){
          storage->t[nt] = MIN(head-BIPARTITE, head2-BIPARTITE);
          storage->h[nt] = MAX(head-BIPARTITE, head2-BIPARTITE);
          storage->w[nt] = WtGETWT(head-BIPARTITE, head2-BIPARTITE, pnwp) + echange;
          nt++;
        }
      }); break;
  default: error("We should never be here.");
  }

  WtChangeStats(nt, storage->t, storage->h, storage->w, pnwp, storage->m);

  memcpy(CHANGE_STAT, storage->m->workspace, N_CHANGE_STATS*sizeof(double));
}

Z_CHANGESTAT_FN(z_on_proj_net){
  GET_AUX_STORAGE(WtNetwork, pnwp);
  GET_STORAGE(StoreWtModelAndWtChanges, storage)

  WtZStats(pnwp, storage->m, FALSE);

  memcpy(CHANGE_STAT, storage->m->workspace, N_CHANGE_STATS*sizeof(double));
}

X_CHANGESTAT_FN(x_on_proj_net){
  GET_AUX_STORAGE(WtNetwork, pnwp);
  GET_STORAGE(StoreWtModelAndWtChanges, storage);
  ModelTerm *_mymtp = mtp;
  WtSEND_X_SIGNAL_INTO(pnwp, storage->m, NULL, _mymtp->dstats, type, data);
}

F_CHANGESTAT_FN(f_on_proj_net){
  GET_AUX_STORAGE(WtNetwork, pnwp);
  GET_STORAGE(StoreWtModelAndWtChanges, storage);

  R_Free(storage->t);
  R_Free(storage->h);
  R_Free(storage->w);
  WtModelDestroy(pnwp, storage->m);
}



I_CHANGESTAT_FN(i__proj_net){
  unsigned int mode = IINPUT_PARAM[0];

  WtNetwork *pnwp = AUX_STORAGE = WtNetworkInitialize(NULL, NULL, NULL, 0, mode == 1 ? BIPARTITE : N_NODES - BIPARTITE, DIRECTED, FALSE);

  switch(mode){
  case 1:
    EXEC_THROUGH_NET_EDGES_PRE(tail, head, e1, {
        EXEC_THROUGH_FINEDGES(head, e2, tail2, {
            if(tail < tail2)
              WtSETWT(tail, tail2, WtGETWT(tail, tail2, pnwp) + 1, pnwp);
          });
      }); break;
  case 2:
    EXEC_THROUGH_NET_EDGES_PRE(tail, head, e1, {
        EXEC_THROUGH_FOUTEDGES(tail, e2, head2, {
            if(head < head2)
              WtSETWT(head-BIPARTITE, head2-BIPARTITE, WtGETWT(head-BIPARTITE, head2-BIPARTITE, pnwp) + 1, pnwp);
          });
      }); break;
  default: error("We should never be here.");
  }
}

U_CHANGESTAT_FN(u__proj_net){
  GET_AUX_STORAGE(WtNetwork, pnwp);
  int echange = edgestate ? -1 : 1;
  unsigned int mode = IINPUT_PARAM[0];

  switch(mode){
  case 1:
    EXEC_THROUGH_FINEDGES(head, e, tail2, {
        if(tail!=tail2)
          WtSETWT(tail, tail2, WtGETWT(tail, tail2, pnwp) + echange, pnwp);
      }); break;
  case 2:
    EXEC_THROUGH_FOUTEDGES(tail, e, head2, {
        if(head!=head2)
          WtSETWT(head-BIPARTITE, head2-BIPARTITE, WtGETWT(head-BIPARTITE, head2-BIPARTITE, pnwp) + echange, pnwp);
      }); break;
  default: error("We should never be here.");
  }
}

F_CHANGESTAT_FN(f__proj_net){
  GET_AUX_STORAGE(WtNetwork, pnwp);
  WtNetworkDestroy(pnwp);
  AUX_STORAGE = NULL;
}
