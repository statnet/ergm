#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_changestat_operator.h"

typedef struct {
  unsigned int nl;
  Network *nwp;
  double *lid;
  double *lmap;
} StoreNetworksAndLayerIDAndLayerMapAndNLayers;

I_CHANGESTAT_FN(i__layer_nets){
  ALLOC_AUX_STORAGE(1, StoreNetworksAndLayerIDAndLayerMapAndNLayers, linfo);
  linfo->nl = INPUT_PARAM[1];
  Vertex lnnodes = N_NODES/linfo->nl, lbip = BIPARTITE/linfo->nl;
  linfo->nwp = Calloc(linfo->nl+1, Network);
  linfo->lid = INPUT_PARAM+2 -1; // The -1 is because Vertex IDs count from 1.
  linfo->lmap = INPUT_PARAM+2+N_NODES -1;
  for(unsigned int l = 1; l <= linfo->nl; l++){
    linfo->nwp[l] = NetworkInitialize(NULL, NULL, 0, lnnodes, DIRECTED, lbip, 0, 0, NULL);
  }
  
  EXEC_THROUGH_NET_EDGES(t, h, e, {
      ToggleEdge(linfo->lmap[t], linfo->lmap[h], linfo->nwp + (int)linfo->lid[t]);
    });
}

U_CHANGESTAT_FN(u__layer_nets){ 
  GET_AUX_STORAGE(StoreNetworksAndLayerIDAndLayerMapAndNLayers, linfo);
  ToggleEdge(linfo->lmap[tail], linfo->lmap[head], linfo->nwp + (int)linfo->lid[tail]);
}

F_CHANGESTAT_FN(f__layer_nets){ 
  GET_AUX_STORAGE(StoreNetworksAndLayerIDAndLayerMapAndNLayers, linfo);
  for(unsigned int l = 1; l <= linfo->nl; l++){
    NetworkDestroy(linfo->nwp + l);
  }
  Free(linfo->nwp);
}

I_CHANGESTAT_FN(i_OnLayer){
  double *inputs = INPUT_PARAM+1; // One for auxiliary storage.
  unsigned int nml = *(inputs++); // Number of layers *in the term*.
  double *lIDmap = inputs; // Which layer does a particular layer in the term refer to?
  inputs += nml;
  
  GET_AUX_STORAGE(StoreNetworksAndLayerIDAndLayerMapAndNLayers, linfo);  
  ALLOC_STORAGE(linfo->nl, Model*, ms);

  for(unsigned int ml=0; ml<nml; ml++){
    ms[ml] = unpack_Model_as_double(&inputs);
    InitStats(linfo->nwp+(int)(lIDmap[ml]), ms[ml]);
  }
}

C_CHANGESTAT_FN(c_OnLayer){
  GET_AUX_STORAGE(StoreNetworksAndLayerIDAndLayerMapAndNLayers, linfo);  
  GET_STORAGE(Model*, ms);
  unsigned int nml = INPUT_PARAM[1];
  double *lIDmap = INPUT_PARAM+2; // Which layer does a particular layer in the model refer to?

  Vertex lt = linfo->lmap[tail], lh = linfo->lmap[head], l = linfo->lid[tail]; 

  // Find the affected model.
  unsigned int ml;
  for(ml=0; ml < nml && lIDmap[ml]!=l; ml++);
  if(ml<nml){ // i.e., affected layer found
    ChangeStats(1, &lt, &lh, linfo->nwp+l, ms[ml]);
    memcpy(CHANGE_STAT, ms[ml]->workspace, N_CHANGE_STATS*sizeof(double));
  }
}

U_CHANGESTAT_FN(u_OnLayer){
  GET_AUX_STORAGE(StoreNetworksAndLayerIDAndLayerMapAndNLayers, linfo);  
  GET_STORAGE(Model*, ms);
  unsigned int nml = INPUT_PARAM[1];
  double *lIDmap = INPUT_PARAM+2; // Which layer does a particular layer in the model refer to?

  Vertex lt = linfo->lmap[tail], lh = linfo->lmap[head], l = linfo->lid[tail]; 
  
  // Find the affected model.
  unsigned int ml;
  for(ml=0; ml < nml && lIDmap[ml]!=l; ml++);
  if(ml<nml){ // i.e., affected layer found
    Model *m = ms[ml];
    UPDATE_STORAGE(lt, lh, linfo->nwp+l, ms[ml], NULL);
  }
}

F_CHANGESTAT_FN(f_OnLayer){
  GET_AUX_STORAGE(StoreNetworksAndLayerIDAndLayerMapAndNLayers, linfo);  
  GET_STORAGE(Model*, ms);
  unsigned int nml = INPUT_PARAM[1];
  double *lIDmap = INPUT_PARAM+2; // Which layer does a particular layer in the model refer to?
  for(unsigned int ml=0; ml<nml; ml++){
    ModelDestroy(linfo->nwp+(int)(lIDmap[ml]), ms[ml]);
  }
}
