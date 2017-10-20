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

  if(DIRECTED){
    ll->symm = inputs - 1; // The -1 is because layer IDs count from 1.
    inputs += ll->nl;
  }else ll->symm = NULL;

  ll->onwp = Calloc(1, Network);
  ll->onwp[0] = NetworkInitialize(NULL, NULL, 0, lnnodes, DIRECTED, lbip, 0, 0, NULL);
  
  /* Set up the layer logic. */

  ll->commands = inputs;
  ll->stacks = Calloc(2*ll->commands[0], double);

  /* Construct the output (logical layer) network: */  
  
  EXEC_THROUGH_NET_EDGES_PRE(t, h, e, {
      Vertex tl = ML_LID_TAIL(ll, t);
      if(ergm_LayerLogic(t, h, ll, 0)){
	ML_SETWT(ll, ML_IO_TAIL(ll, t), ML_IO_HEAD(ll, h), 1);
      }
      // If the physical layer is symmetrized then also check the
      // layer logic of its recirpocation.
      if(ll->symm && ll->symm[tl]!=0 && ergm_LayerLogic(h, t, ll, 0)){
	ML_SETWT(ll, ML_IO_HEAD(ll, h), ML_IO_TAIL(ll, t), 1);	
      }
    });
}

U_CHANGESTAT_FN(u__layer_net){ 
  GET_AUX_STORAGE(StoreLayerLogic, ll);
  Vertex tl = ML_LID_TAIL(ll, tail);
  if(ergm_LayerLogic(tail, head, ll, TRUE)){
    ML_TOGGLE(ll, ML_IO_TAIL(ll, tail), ML_IO_HEAD(ll, head));
  }
  if(ll->symm && ll->symm[tl]!=0 && ergm_LayerLogic(head, tail, ll, TRUE)){
    ML_TOGGLE(ll, ML_IO_HEAD(ll, head), ML_IO_TAIL(ll, tail));
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
  
/*   EXEC_THROUGH_NET_EDGES_PRE(t, h, e, { */
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
    double *inputs = INPUT_ATTRIB+nml+1; // Rewind to the start of model spec.
    ms[ml] = unpack_Model_as_double(&inputs);
    InitStats(ll->onwp, ms[ml]);
  }
}

C_CHANGESTAT_FN(c_OnLayer){
  GET_STORAGE(Model*, ms);
  unsigned int nml = *INPUT_ATTRIB;
  double *w = INPUT_ATTRIB+1;

  // Find the affected models.
  for(unsigned int ml=0; ml < nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex tl = ML_LID_TAIL(ll, tail);
    if(!ll->symm || ll->symm[tl]==0){ // Toggle not symmetrized
      if(ergm_LayerLogic(tail, head, ll, TRUE)){ // network affected
	Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
	ChangeStats(1, &lt, &lh, ll->onwp, ms[ml]);
	for(unsigned int i=0; i<N_CHANGE_STATS; i++)
	  CHANGE_STAT[i] += ms[ml]->workspace[i] * w[ml];
      }
    }else{ // Toggle symmetrized
      unsigned int nt = 0;
      Vertex lt[2], lh[2];
      if(ergm_LayerLogic(tail, head, ll, TRUE)){
	lt[nt] = ML_IO_TAIL(ll, tail);
	lh[nt] = ML_IO_HEAD(ll, head);
	nt++;
      }
      if(ergm_LayerLogic(head, tail, ll, TRUE)){
	// FIXME: This might break bipartite networks. We do not have
	// directed bipartite networks at this time, but this may
	// change in the future.
	lh[nt] = ML_IO_TAIL(ll, tail);
	lt[nt] = ML_IO_HEAD(ll, head);
	nt++;
      }
      if(nt){
	ChangeStats(nt, lt, lh, ll->onwp, ms[ml]);
	for(unsigned int i=0; i<N_CHANGE_STATS; i++)
	  CHANGE_STAT[i] += ms[ml]->workspace[i] * w[ml];
      }
    }
  }
}

U_CHANGESTAT_FN(u_OnLayer){
  GET_STORAGE(Model*, ms);
  unsigned int nml = *INPUT_ATTRIB;

  // Find the affected models.
  for(unsigned int ml=0; ml < nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex tl = ML_LID_TAIL(ll, tail);
    if(!ll->symm || ll->symm[tl]==0){ // Toggle not symmetrized
      if(ergm_LayerLogic(tail, head, ll, TRUE)){ // network affected
	Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
	UPDATE_STORAGE(lt, lh, ll->onwp, ms[ml], NULL);
      }
    }else{ // Toggle symmetrized
      Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
      if(ergm_LayerLogic(tail, head, ll, TRUE)){
	UPDATE_STORAGE(lt, lh, ll->onwp, ms[ml], NULL);
      }	
      if(ergm_LayerLogic(head, tail, ll, TRUE)){
	// We need to make a provisional toggle here, since
	// UPDATE_STORAGE expects one toggle at a time.
	ToggleEdge(lt, lh, ll->onwp);
	// FIXME: This might break bipartite networks. We do not have
	// directed bipartite networks at this time, but this may
	// change in the future.
	Vertex lt = ML_IO_TAIL(ll, head), lh = ML_IO_HEAD(ll, tail);
	UPDATE_STORAGE(lt, lh, ll->onwp, ms[ml], NULL);
	// Reverse the provisional toggle.
	ToggleEdge(lt, lh, ll->onwp);
      }
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



/*****************
 changestat: d_ldegree_ML_sum
*****************/
C_CHANGESTAT_FN(c_ldegree_ML_sum) { 
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double *dirs = inputs;
  inputs += nml;
  double *degs = inputs;

  /* *** don't forget tail -> head */

  unsigned int taildeg = 0, headdeg = 0;
  int tdegchange = 0, hdegchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *id=ML_IN_DEG(ll), *od=ML_OUT_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    int degchange = ergm_LayerLogic(tail, head, ll, TRUE);
    if(dirs[ml]==0 || dirs[ml]==+1){
      tdegchange += degchange;
      taildeg += od[lt];
    }
    if(dirs[ml]==0 || dirs[ml]==-1){
      hdegchange += degchange;
      headdeg += id[lh];
    }
  }
  
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
    Vertex deg = (Vertex)degs[j];
    CHANGE_STAT[j] += (taildeg + tdegchange == deg) - (taildeg == deg);
    CHANGE_STAT[j] += (headdeg + hdegchange == deg) - (headdeg == deg);
  }
}

/*****************
 changestat: d_ldegree_by_attr_ML_sum
*****************/
C_CHANGESTAT_FN(c_ldegree_by_attr_ML_sum) { 
  /* The inputparams are assumed to be set up as follows:
  The first 2*nstats values are in pairs:  (degree, attrvalue)
  The values following the first 2*nstats values are the nodal attributes.
  */
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double *dirs = inputs;
  inputs += nml;

  /* *** don't forget tail -> head */    

  unsigned int taildeg = 0, headdeg = 0;
  int tdegchange = 0, hdegchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *id=ML_IN_DEG(ll), *od=ML_OUT_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    int degchange = ergm_LayerLogic(tail, head, ll, TRUE);
    if(dirs[ml]==0 || dirs[ml]==+1){
      tdegchange += degchange;
      taildeg += od[lt];
    }
    if(dirs[ml]==0 || dirs[ml]==-1){
      hdegchange += degchange;
      headdeg += id[lh];
    }
  }   

  int tailattr = inputs[2*N_CHANGE_STATS + tail - 1]; 
  int headattr = inputs[2*N_CHANGE_STATS + head - 1]; 
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
    Vertex d = (Vertex)inputs[2*j];
    int testattr = inputs[2*j + 1]; 
    if (tailattr == testattr)  /* we have tail attr match */
      CHANGE_STAT[j] += (taildeg + tdegchange == d) - (taildeg == d);
    if (headattr == testattr)  /* we have head attr match */
      CHANGE_STAT[j] += (headdeg + hdegchange == d) - (headdeg == d);
  }
}

/*****************
 changestat: d_ldegree_ML_sum
*****************/
C_CHANGESTAT_FN(c_gwldegree_ML_sum) { 
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double *dirs = inputs;
  inputs += nml;
  double decay = *(inputs++);
  double oneexpd = 1.0-exp(-decay);

  /* *** don't forget tail -> head */

  unsigned int taildeg = 0, headdeg = 0;
  int tdegchange = 0, hdegchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *id=ML_IN_DEG(ll), *od=ML_OUT_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    int degchange = ergm_LayerLogic(tail, head, ll, TRUE);
    if(dirs[ml]==0 || dirs[ml]==+1){
      tdegchange += degchange;
      taildeg += od[lt];
    }
    if(dirs[ml]==0 || dirs[ml]==-1){
      hdegchange += degchange;
      headdeg += id[lh];
    }
  }
  
  CHANGE_STAT[0] =
    exp(decay) * (
		  ((1-pow(oneexpd,taildeg+tdegchange)) - (1-pow(oneexpd,taildeg))) +
		  ((1-pow(oneexpd,headdeg+hdegchange)) - (1-pow(oneexpd,headdeg)))
		  );
}

/*****************
 changestat: d_ldegree_by_attr_ML_sum
*****************/
C_CHANGESTAT_FN(c_gwldegree_by_attr_ML_sum) { 
  /* The inputparams are assumed to be set up as follows:
  The first 2*nstats values are in pairs:  (degree, attrvalue)
  The values following the first 2*nstats values are the nodal attributes.
  */
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double *dirs = inputs;
  inputs += nml;
  double decay = *(inputs++);
  double oneexpd = 1.0-exp(-decay);

  /* *** don't forget tail -> head */    

  unsigned int taildeg = 0, headdeg = 0;
  int tdegchange = 0, hdegchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *id=ML_IN_DEG(ll), *od=ML_OUT_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    int degchange = ergm_LayerLogic(tail, head, ll, TRUE);
    if(dirs[ml]==0 || dirs[ml]==+1){
      tdegchange += degchange;
      taildeg += od[lt];
    }
    if(dirs[ml]==0 || dirs[ml]==-1){
      hdegchange += degchange;
      headdeg += id[lh];
    }
  }   

  int tailattr = inputs[2*N_CHANGE_STATS + tail - 1]; 
  int headattr = inputs[2*N_CHANGE_STATS + head - 1]; 
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
    int testattr = inputs[2*j + 1]; 
    if (tailattr == testattr)  /* we have tail attr match */
      CHANGE_STAT[j] +=
	exp(decay) * ((1-pow(oneexpd,taildeg+tdegchange)) - (1-pow(oneexpd,taildeg)));
    if (headattr == testattr)  /* we have head attr match */
      CHANGE_STAT[j] +=
	exp(decay) * ((1-pow(oneexpd,headdeg+hdegchange)) - (1-pow(oneexpd,headdeg)));

  }
}
