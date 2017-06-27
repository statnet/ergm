#include "ergm_changestat_multilayer.h"
#include "ergm_changestat.h"
#include "ergm_changestat_operator.h"


I_CHANGESTAT_FN(i__layer_net){
  double *inputs = INPUT_PARAM;
  ALLOC_AUX_STORAGE(1, StoreLayerLogic, ll); inputs++;
  ll->nl = *(inputs++);
  ll->inwp = nwp;

  /* Set up the layer information. */
  ll->lid = inputs - 1; // The -1 is because Vertex IDs count from 1.
  inputs += N_NODES;
  ll->lmap = inputs - 1;
  inputs += N_NODES;

  Vertex lnnodes, lbip;
  if(BIPARTITE){
    lbip = lnnodes = *(inputs++);
    lnnodes += *(inputs++);
    inputs += (ll->nl-1)*2; // There will be a total of nl*2 network sizes.
  }else{
    lbip = 0;
    lnnodes = *(inputs++);
    inputs += (ll->nl-1); // There will be a total of nl network sizes.
  }

  ll->onwp = Calloc(1, Network);
  ll->onwp[0] = NetworkInitialize(NULL, NULL, 0, lnnodes, DIRECTED, lbip, 0, 0, NULL);
  
  /* Set up the layer logic. */

  ll->commands = inputs;
  ll->stacks = Calloc(2*ll->commands[0], double);

  /* Construct the output (logical layer) network: */  
  
  EXEC_THROUGH_NET_EDGES(t, h, e, {
      if(ergm_LayerLogic(t, h, ll, TRUE)){
	ML_TOGGLE(ll, ML_IO_TAIL(ll, t), ML_IO_HEAD(ll, h));
      }
    });
}

U_CHANGESTAT_FN(u__layer_net){ 
  GET_AUX_STORAGE(StoreLayerLogic, ll);
  if(ergm_LayerLogic(tail, head, ll, TRUE)){
    ML_TOGGLE(ll, ML_IO_TAIL(ll, tail), ML_IO_HEAD(ll, head));
  }
}

F_CHANGESTAT_FN(f__layer_net){ 
  GET_AUX_STORAGE(StoreLayerLogic, ll);
  NetworkDestroy(ll->onwp);
  Free(ll->onwp);
  Free(ll->stacks);
}

/* I_CHANGESTAT_FN(i__layer_nets){ */
/*   ALLOC_AUX_STORAGE(1, StoreNetsAndLIDAndLMapAndNL, li); */
/*   li->nl = INPUT_PARAM[1]; */
/*   Vertex lnnodes = N_NODES/li->nl, lbip = BIPARTITE/li->nl; */
/*   li->nwp = Calloc(li->nl+1, Network); */
/*   li->lid = INPUT_PARAM+2 -1; // The -1 is because Vertex IDs count from 1. */
/*   li->lmap = INPUT_PARAM+2+N_NODES -1; */
/*   for(unsigned int l = 1; l <= li->nl; l++){ */
/*     li->nwp[l] = NetworkInitialize(NULL, NULL, 0, lnnodes, DIRECTED, lbip, 0, 0, NULL); */
/*   } */
  
/*   EXEC_THROUGH_NET_EDGES(t, h, e, { */
/*       ToggleEdge(li->lmap[t], li->lmap[h], li->nwp + (int)li->lid[t]); */
/*     }); */
/* } */

/* U_CHANGESTAT_FN(u__layer_nets){  */
/*   GET_AUX_STORAGE(StoreNetsAndLIDAndLMapAndNL, li); */
/*   ToggleEdge(li->lmap[tail], li->lmap[head], li->nwp + (int)li->lid[tail]); */
/* } */

/* F_CHANGESTAT_FN(f__layer_nets){  */
/*   GET_AUX_STORAGE(StoreNetsAndLIDAndLMapAndNL, li); */
/*   for(unsigned int l = 1; l <= li->nl; l++){ */
/*     NetworkDestroy(li->nwp + l); */
/*   } */
/*   Free(li->nwp); */
/* } */

I_CHANGESTAT_FN(i_OnLayer){
  
  unsigned int nml = *INPUT_ATTRIB; // Number of layers *in the term*. Already shifted past the auxiliaries.
  
  ALLOC_STORAGE(nml, Model*, ms);

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    double *inputs = INPUT_ATTRIB+1; // Rewind to the start of model spec.
    ms[ml] = unpack_Model_as_double(&inputs);
    InitStats(ll->onwp, ms[ml]);
  }
}

C_CHANGESTAT_FN(c_OnLayer){
  GET_STORAGE(Model*, ms);
  unsigned int nml = *INPUT_ATTRIB;

  // Find the affected models.
  for(unsigned int ml=0; ml < nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    if(ergm_LayerLogic(tail, head, ll, TRUE)){ // network affected
      Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
      ChangeStats(1, &lt, &lh, ll->onwp, ms[ml]);
      for(unsigned int i=0; i<N_CHANGE_STATS; i++)
	CHANGE_STAT[i] += ms[ml]->workspace[i];
    }
  }
}

U_CHANGESTAT_FN(u_OnLayer){
  GET_STORAGE(Model*, ms);
  unsigned int nml = *INPUT_ATTRIB;

  // Find the affected models.
  for(unsigned int ml=0; ml < nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    if(ergm_LayerLogic(tail, head, ll, TRUE)){ // network affected
      Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
      Model *m = ms[ml];
      UPDATE_STORAGE(lt, lh, ll->onwp, ms[ml], NULL);
    }
  }
}

F_CHANGESTAT_FN(f_OnLayer){
  GET_STORAGE(Model*, ms);
  unsigned int nml = *INPUT_ATTRIB;
  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    ModelDestroy(ll->onwp, ms[ml]);
  }
}

/* layerCMB: Conway-Maxwell-Binomial for the sum of layer combinations */

C_CHANGESTAT_FN(c_layerCMB){
  unsigned int nml = *INPUT_ATTRIB;

  // FIXME: Cache current values, perhaps via a valued auxiliary?

  unsigned int oldct=0, newct=0;
  for(unsigned int ml=0; ml < nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    unsigned int v = ergm_LayerLogic(tail, head, ll, 2);
    if(v&1) oldct++; // Pre-toggle edge present.
    if(v&2) newct++; // Post-toggle edge present.
  }
  
  CHANGE_STAT[0] = lgamma1p(newct)-lgamma1p(oldct) + lgamma1p(nml-newct)-lgamma1p(nml-oldct); 
}
