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
    unsigned int need_symm = FALSE;
    for(unsigned int l=1; l<=ll->nl; l++){
      if(ll->symm[l]){
	need_symm = TRUE;
	break;
      }
    }
    if(!need_symm) ll->symm = NULL;
    
    inputs += ll->nl;
  }else ll->symm = NULL;

  ll->onwp = Calloc(1, Network);
  ll->onwp[0] = NetworkInitialize(NULL, NULL, 0, lnnodes, DIRECTED, lbip, 0, 0, NULL);
  
  /* Set up the layer logic. */

  ll->commands = inputs;
  ll->stacks = Calloc(2*ll->commands[0], double);

  /* Construct the output (logical layer) network: */  

  EXEC_THROUGH_NET_EDGES_PRE(t, h, e, {
      Vertex at[2];
      Vertex ah[2];
      unsigned int nt = ergm_LayerLogic_affects(t, h, ll, 0, at, ah);
      for(unsigned int i=0; i<nt; i++){
	ML_SETWT(ll, at[i], ah[i], 1);
      }
    });
}

U_CHANGESTAT_FN(u__layer_net){ 
  GET_AUX_STORAGE(StoreLayerLogic, ll);
  Vertex at[2], ah[2];
  unsigned int nt = ergm_LayerLogic_affects(tail, head, ll, 1, at, ah);
  for(unsigned int i=0; i<nt; i++){
    ML_TOGGLE(ll, at[i], ah[i]);
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
    Vertex at[2], ah[2];
    unsigned int nt = ergm_LayerLogic_affects(tail, head, ll, 1, at, ah);
    if(nt){
      ChangeStats(nt, at, ah, ll->onwp, ms[ml]);
      for(unsigned int i=0; i<N_CHANGE_STATS; i++)
	CHANGE_STAT[i] += ms[ml]->workspace[i] * w[ml];
    }
  }
}

U_CHANGESTAT_FN(u_OnLayer){
  GET_STORAGE(Model*, ms);
  unsigned int nml = *INPUT_ATTRIB;

  // Find the affected models.
  for(unsigned int ml=0; ml < nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex at[2], ah[2];
    unsigned int nt = ergm_LayerLogic_affects(tail, head, ll, 1, at, ah);
    int i;
    for(i=0; i<nt; i++){
      UPDATE_STORAGE(at[i], ah[i], ll->onwp, ms[ml], NULL);
      // We need to make a provisional toggle here, since
      // UPDATE_STORAGE expects one toggle at a time. Note that
      // ll->onwp is "owned" by the .layer.net auxiliary, so this may
      // cause a race condition for multithreaded term evaluation. The
      // "correct" solution might involve duplicating the network,
      // though some locking mechanism might also work.
      if(i+1 < nt) ToggleEdge(at[i], ah[i], ll->onwp);
    }
    // Reverse the provisional toggle.
    i--; while(--i>=0) ToggleEdge(at[i], ah[i], ll->onwp);
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

  unsigned int oldct_th=0, newct_th=0,
    oldct_ht=0, newct_ht=0;
  for(unsigned int ml=0; ml < nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    unsigned int v = ergm_LayerLogic2(lt, lh, tail, head, ll, 2);
    if(v&1) oldct_th++; // Pre-toggle edge present.
    if(v&2) newct_th++; // Post-toggle edge present.
    
    v = ergm_LayerLogic2(lh, lt, tail, head, ll, 2);
    if(v&1) oldct_ht++; // Pre-toggle edge present.
    if(v&2) newct_ht++; // Post-toggle edge present.
  }
  
  CHANGE_STAT[0] =
    +(newct_th!=oldct_th? lgamma1p(newct_th)-lgamma1p(oldct_th) + lgamma1p(nml-newct_th)-lgamma1p(nml-oldct_th) : 0)
    +(newct_ht!=oldct_ht? lgamma1p(newct_ht)-lgamma1p(oldct_ht) + lgamma1p(nml-newct_ht)-lgamma1p(nml-oldct_ht) : 0);    
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
