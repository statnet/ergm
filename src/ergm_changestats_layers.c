#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_storage.h"

typedef struct {
  Network *nwp;
  double *lid;
  double *lmap;
} StoreNetworksAndLayerIDAndLayerMap;

I_CHANGESTAT_FN(i__layer_nets){
  ALLOC_AUX_STORAGE(1, StoreNetworksAndLayerIDAndLayerMap, linfo);
  unsigned int nl = INPUT_PARAM[0];
  Vertex lnnodes = N_NODES/nl, lbip = BIPARTITE/nl;
  linfo->nwp = Calloc(nl+1, Network);
  linfo->lid = INPUT_PARAM+1;
  linfo->lmap = INPUT_PARAM+1+N_NODES;
  for(unsigned int l = 1; l <= nl; l++){
    linfo->nwp[l] = NetworkInitialize(NULL, NULL, 0, lnnodes, DIRECTED, lbip, 0, 0, NULL);
  }
  
  EXEC_THROUGH_NET_EDGES(t, h, e, {
      ToggleEdge(linfo->lmap[t], linfo->lmap[h], linfo->nwp + (int)linfo->lid[t]);
    });
}

U_CHANGESTAT_FN(u__layer_nets){ 
  GET_AUX_STORAGE(StoreNetworksAndLayerIDAndLayerMap, linfo);
  ToggleEdge(linfo->lmap[tail], linfo->lmap[head], linfo->nwp + (int)linfo->lid[tail]);
}

F_CHANGESTAT_FN(f__layer_nets){ 
  unsigned int nl = INPUT_PARAM[0];
  GET_AUX_STORAGE(StoreNetworksAndLayerIDAndLayerMap, linfo);
  for(unsigned int l = 1; l <= nl; l++){
    NetworkDestroy(linfo->nwp + l);
  }
}
