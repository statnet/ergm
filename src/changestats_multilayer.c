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

  ll->onwp = NetworkInitialize(NULL, NULL, 0, lnnodes, DIRECTED, lbip, 0, 0, NULL);
  
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

  unsigned int todeg = 0, hodeg = 0, tideg = 0, hideg = 0;
  int todegchange = 0, hodegchange = 0, tidegchange = 0, hidegchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *id=ML_IN_DEG(ll), *od=ML_OUT_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    int degchange_th = ergm_LayerLogic2(lt, lh, tail, head, ll, TRUE),
      degchange_ht = ergm_LayerLogic2(lh, lt, tail, head, ll, TRUE);
    if(dirs[ml]==0 || dirs[ml]==+1){
      todegchange += degchange_th;
      hodegchange += degchange_ht;
      todeg += od[lt];
      hodeg += od[lh];
    }
    if(dirs[ml]==0 || dirs[ml]==-1){
      hidegchange += degchange_th;
      tidegchange += degchange_ht;
      hideg += id[lh];
      tideg += id[lt];
    }
  }
  
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
    Vertex deg = (Vertex)degs[j];
    CHANGE_STAT[j] += (todeg + todegchange == deg) - (todeg == deg)
      + (tideg + tidegchange == deg) - (tideg == deg) 
      + (hodeg + hodegchange == deg) - (hodeg == deg)
      + (hideg + hidegchange == deg) - (hideg == deg);
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

  unsigned int todeg = 0, hodeg = 0, tideg = 0, hideg = 0;
  int todegchange = 0, hodegchange = 0, tidegchange = 0, hidegchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *id=ML_IN_DEG(ll), *od=ML_OUT_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    int degchange_th = ergm_LayerLogic2(lt, lh, tail, head, ll, TRUE),
      degchange_ht = ergm_LayerLogic2(lh, lt, tail, head, ll, TRUE);
    if(dirs[ml]==0 || dirs[ml]==+1){
      todegchange += degchange_th;
      hodegchange += degchange_ht;
      todeg += od[lt];
      hodeg += od[lh];
    }
    if(dirs[ml]==0 || dirs[ml]==-1){
      hidegchange += degchange_th;
      tidegchange += degchange_ht;
      hideg += id[lh];
      tideg += id[lt];
    }
  }   

  int tailattr = inputs[2*N_CHANGE_STATS + tail - 1]; 
  int headattr = inputs[2*N_CHANGE_STATS + head - 1]; 
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
    Vertex deg = (Vertex)inputs[2*j];
    int testattr = inputs[2*j + 1]; 
    if (tailattr == testattr)  /* we have tail attr match */
      CHANGE_STAT[j] += (todeg + todegchange == deg) - (todeg == deg)
	+ (tideg + tidegchange == deg) - (tideg == deg);
    if (headattr == testattr)  /* we have head attr match */
      CHANGE_STAT[j] += (hodeg + hodegchange == deg) - (hodeg == deg)
	+ (hideg + hidegchange == deg) - (hideg == deg);
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

  unsigned int todeg = 0, hodeg = 0, tideg = 0, hideg = 0;
  int todegchange = 0, hodegchange = 0, tidegchange = 0, hidegchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *id=ML_IN_DEG(ll), *od=ML_OUT_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    int degchange_th = ergm_LayerLogic2(lt, lh, tail, head, ll, TRUE),
      degchange_ht = ergm_LayerLogic2(lh, lt, tail, head, ll, TRUE);
    if(dirs[ml]==0 || dirs[ml]==+1){
      todegchange += degchange_th;
      hodegchange += degchange_ht;
      todeg += od[lt];
      hodeg += od[lh];
    }
    if(dirs[ml]==0 || dirs[ml]==-1){
      hidegchange += degchange_th;
      tidegchange += degchange_ht;
      hideg += id[lh];
      tideg += id[lt];
    }
  }
  
  CHANGE_STAT[0] =
    exp(decay) * (
		  ((1-pow(oneexpd,todeg+todegchange)) - (1-pow(oneexpd,todeg))) +
		  ((1-pow(oneexpd,hodeg+hodegchange)) - (1-pow(oneexpd,hodeg))) +
		  ((1-pow(oneexpd,tideg+tidegchange)) - (1-pow(oneexpd,tideg))) +
		  ((1-pow(oneexpd,hideg+hidegchange)) - (1-pow(oneexpd,hideg)))
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

  unsigned int todeg = 0, hodeg = 0, tideg = 0, hideg = 0;
  int todegchange = 0, hodegchange = 0, tidegchange = 0, hidegchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *id=ML_IN_DEG(ll), *od=ML_OUT_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    int degchange_th = ergm_LayerLogic2(lt, lh, tail, head, ll, TRUE),
      degchange_ht = ergm_LayerLogic2(lh, lt, tail, head, ll, TRUE);
    if(dirs[ml]==0 || dirs[ml]==+1){
      todegchange += degchange_th;
      hodegchange += degchange_ht;
      todeg += od[lt];
      hodeg += od[lh];
    }
    if(dirs[ml]==0 || dirs[ml]==-1){
      hidegchange += degchange_th;
      tidegchange += degchange_ht;
      hideg += id[lh];
      tideg += id[lt];
    }
  }   

  int tailattr = inputs[2*N_CHANGE_STATS + tail - 1]; 
  int headattr = inputs[2*N_CHANGE_STATS + head - 1]; 
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
    int testattr = inputs[2*j + 1]; 
    if (tailattr == testattr)  /* we have tail attr match */
      CHANGE_STAT[j] +=
	    exp(decay) * (
			  ((1-pow(oneexpd,todeg+todegchange)) - (1-pow(oneexpd,todeg))) +
			  ((1-pow(oneexpd,tideg+tidegchange)) - (1-pow(oneexpd,tideg)))
			  );
    if (headattr == testattr)  /* we have head attr match */
      CHANGE_STAT[j] +=
	    exp(decay) * (
			  ((1-pow(oneexpd,hodeg+hodegchange)) - (1-pow(oneexpd,hodeg))) +
			  ((1-pow(oneexpd,hideg+hidegchange)) - (1-pow(oneexpd,hideg)))
			  );
  }
}

/*****************
 changestat: c_twopathL
*****************/
C_CHANGESTAT_FN(c_twostarL) { 
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int typeID = inputs[0];
  unsigned int distinct = inputs[1];
  
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll1, 0);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll2, 1);
  
  Vertex lt = ML_IO_TAIL(ll1, tail), lh = ML_IO_HEAD(ll1, head);

  int
    change1_th = ergm_LayerLogic2(lt, lh, tail, head, ll1, TRUE),
    change1_ht = ergm_LayerLogic2(lh, lt, tail, head, ll1, TRUE),
    change2_th = ergm_LayerLogic2(lt, lh, tail, head, ll2, TRUE),
    change2_ht = ergm_LayerLogic2(lh, lt, tail, head, ll2, TRUE);


  // Need int here since we need signed arithmetic.
  int *od1 = (int*) ML_OUT_DEG(ll1), *od2 = (int*) ML_OUT_DEG(ll2),
    *id1 = (int*) ML_IN_DEG(ll1), *id2 = (int*) ML_IN_DEG(ll2);
  
  switch(typeID){
  case 1: // out
    if(change1_th || change2_th){ // lt's outstar counts change.
      
      // Calculate the change in number of relevant coincident relations.
      int change12 = 0;
      if(distinct){
	if(ML_IS_OUTEDGE(ll1, lt, lh)) change12 += change2_th; 
	if(ML_IS_OUTEDGE(ll2, lt, lh)) change12 += change1_th;
	change12 += change1_th * change2_th;
      }
	
      CHANGE_STAT[0] += (od1[lt]+change1_th)*(od2[lt]+change2_th) - od1[lt]*od2[lt] - change12;
    }
    
    if(change1_ht || change2_ht){ // lh's outstar counts change.

      // Calculate the change in number of relevant coincident relations.
      int change12 = 0;
      if(distinct){
	if(ML_IS_OUTEDGE(ll1, lh, lt)) change12 += change2_ht; 
	if(ML_IS_OUTEDGE(ll2, lh, lt)) change12 += change1_ht;
	change12 += change1_ht * change2_ht;
      }

      CHANGE_STAT[0] += (od1[lh]+change1_ht)*(od2[lh]+change2_ht) - od1[lh]*od2[lh] - change12;
    }
    break;
  case 2: // in
    if(change1_th || change2_th){ // lh's instar counts change.
      
      // Calculate the change in number of relevant coincident relations.
      int change12 = 0;
      if(distinct){
	if(ML_IS_OUTEDGE(ll1, lt, lh)) change12 += change2_th; 
	if(ML_IS_OUTEDGE(ll2, lt, lh)) change12 += change1_th;
	change12 += change1_th * change2_th;
      }

      CHANGE_STAT[0] += (id1[lh]+change1_th)*(id2[lh]+change2_th) - id1[lh]*id2[lh] - change12;
    }
    
    if(change1_ht || change2_ht){ // lt's instar counts change.

      // Calculate the change in number of relevant coincident relations.
      int change12 = 0;
      if(distinct){
	if(ML_IS_OUTEDGE(ll1, lh, lt)) change12 += change2_ht; 
	if(ML_IS_OUTEDGE(ll2, lh, lt)) change12 += change1_ht;
	change12 += change1_ht * change2_ht;
      }

      CHANGE_STAT[0] += (id1[lt]+change1_ht)*(id2[lt]+change2_ht) - id1[lt]*id2[lt] - change12;
    }
    break;
  case 3: // path
    if(change1_ht || change2_th){ // two-path count through lt changes.

      // Calculate the change in number of relevant coincident relations.
      int change12 = 0;
      if(distinct){
	if(ML_IS_OUTEDGE(ll1, lh, lt)) change12 += change2_th; 
	if(ML_IS_OUTEDGE(ll2, lt, lh)) change12 += change1_ht;
	change12 += change1_th * change2_ht;
      }

      CHANGE_STAT[0] += (id1[lt]+change1_ht)*(od2[lt]+change2_th) - id1[lt]*od2[lt] - change12;
    }
    
    if(change1_th || change2_ht){ // two-path count through lh changes.

      // Calculate the change in number of relevant coincident relations.
      int change12 = 0;
      if(distinct){
	if(ML_IS_INEDGE(ll1, lh, lt)) change12 += change2_ht; 
	if(ML_IS_OUTEDGE(ll2, lh, lt)) change12 += change1_th;
	change12 += change1_th * change2_ht;
      }

      CHANGE_STAT[0] += (id1[lh]+change1_th)*(od2[lh]+change2_ht) - id1[lh]*od2[lh] - change12;
    }
    break;
  }
}

/*****************
 changestat: d_b1degree_ML_sum
*****************/
C_CHANGESTAT_FN(c_b1degree_ML_sum) { 
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double *degs = inputs;

  /* *** don't forget tail -> head */

  unsigned int b1deg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *od=ML_OUT_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    b1deg += od[lt];
  }
  
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
    Vertex deg = (Vertex)degs[j];
    CHANGE_STAT[j] += (b1deg + degchange == deg) - (b1deg == deg);
  }
}

/*****************
 changestat: d_b1degree_by_attr_ML_sum
*****************/
C_CHANGESTAT_FN(c_b1degree_by_attr_ML_sum) { 
  /* The inputparams are assumed to be set up as follows:
  The first 2*nstats values are in pairs:  (degree, attrvalue)
  The values following the first 2*nstats values are the nodal attributes.
  */
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);

  /* *** don't forget tail -> head */    

  unsigned int b1deg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *od=ML_IN_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    b1deg += od[lt];
  }   

  int tailattr = inputs[2*N_CHANGE_STATS + head - 1]; 
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
    Vertex d = (Vertex)inputs[2*j];
    int testattr = inputs[2*j + 1]; 
    if (tailattr == testattr)  /* we have head attr match */
      CHANGE_STAT[j] += (b1deg + degchange == d) - (b1deg == d);
  }
}

/*****************
 changestat: d_b2degree_ML_sum
*****************/
C_CHANGESTAT_FN(c_b2degree_ML_sum) { 
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double *degs = inputs;

  /* *** don't forget tail -> head */

  unsigned int b2deg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *id=ML_IN_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    b2deg += id[lh];
  }
  
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
    Vertex deg = (Vertex)degs[j];
    CHANGE_STAT[j] += (b2deg + degchange == deg) - (b2deg == deg);
  }
}

/*****************
 changestat: d_b2degree_by_attr_ML_sum
*****************/
C_CHANGESTAT_FN(c_b2degree_by_attr_ML_sum) { 
  /* The inputparams are assumed to be set up as follows:
  The first 2*nstats values are in pairs:  (degree, attrvalue)
  The values following the first 2*nstats values are the nodal attributes.
  */
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);

  /* *** don't forget tail -> head */    

  unsigned int b2deg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *id=ML_IN_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    b2deg += id[lh];
  }   

  int headattr = inputs[2*N_CHANGE_STATS + head - 1]; 
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
    Vertex d = (Vertex)inputs[2*j];
    int testattr = inputs[2*j + 1]; 
    if (headattr == testattr)  /* we have head attr match */
      CHANGE_STAT[j] += (b2deg + degchange == d) - (b2deg == d);
  }
}

/*****************
 changestat: d_degree_ML_sum
*****************/
C_CHANGESTAT_FN(c_degree_ML_sum) { 
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double *degs = inputs;

  /* *** don't forget tail -> head */

  unsigned int taildeg = 0, headdeg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *id=ML_IN_DEG(ll), *od=ML_OUT_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    taildeg += od[lt] + id[lt];
    headdeg += od[lh] + id[lh];
  }
  
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
    Vertex deg = (Vertex)degs[j];
    CHANGE_STAT[j] += (taildeg + degchange == deg) - (taildeg == deg);
    CHANGE_STAT[j] += (headdeg + degchange == deg) - (headdeg == deg);
  }
}

/*****************
 changestat: d_degree_by_attr_ML_sum
*****************/
C_CHANGESTAT_FN(c_degree_by_attr_ML_sum) { 
  /* The inputparams are assumed to be set up as follows:
  The first 2*nstats values are in pairs:  (degree, attrvalue)
  The values following the first 2*nstats values are the nodal attributes.
  */
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);

  /* *** don't forget tail -> head */    

  unsigned int taildeg = 0, headdeg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *id=ML_IN_DEG(ll), *od=ML_OUT_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    taildeg += od[lt] + id[lt];
    headdeg += od[lh] + id[lh];
  }   

  int tailattr = inputs[2*N_CHANGE_STATS + tail - 1]; 
  int headattr = inputs[2*N_CHANGE_STATS + head - 1]; 
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
    Vertex d = (Vertex)inputs[2*j];
    int testattr = inputs[2*j + 1]; 
    if (tailattr == testattr)  /* we have tail attr match */
      CHANGE_STAT[j] += (taildeg + degchange == d) - (taildeg == d);
    if (headattr == testattr)  /* we have head attr match */
      CHANGE_STAT[j] += (headdeg + degchange == d) - (headdeg == d);
  }
}

/*****************
 changestat: d_gwb1degree_ML_sum
*****************/
C_CHANGESTAT_FN(c_gwb1degree_ML_sum) { 
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double decay = *(inputs++);
  double oneexpd = 1.0-exp(-decay);

  unsigned int b1deg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *od=ML_OUT_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    b1deg += od[lt];
  }

    
  /* *** don't forget tail -> head */    
  CHANGE_STAT[0] =
    exp(decay) * ((1-pow(oneexpd,b1deg+degchange)) - (1-pow(oneexpd,b1deg)));
}

/*****************
 changestat: d_gwb1degree_by_attr_ML_sum
*****************/
C_CHANGESTAT_FN(c_gwb1degree_by_attr_ML_sum) { 
  /*The inputparams are assumed to be set up as follows:
    The first value is the decay parameter (as in Hunter et al, JASA 200?)
    The next sequence of values is the nodal attributes, coded as integers
         from 1 through N_CHANGE_STATS
  */
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double decay = *(inputs++);
  double *attrs = inputs-1;
  double oneexpd = 1.0-exp(-decay);

  unsigned int b1deg = 0, headdeg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *od=ML_OUT_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    b1deg += od[lt];
  }

    
  /* *** don't forget tail -> head */    
  CHANGE_STAT[(int)attrs[tail]] =
    exp(decay) * ((1-pow(oneexpd,b1deg+degchange)) - (1-pow(oneexpd,b1deg)));
}

/*****************
 changestat: d_gwdegree_ML
*****************/
C_CHANGESTAT_FN(c_gwdegree_ML_sum) { 
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double decay = *(inputs++);
  double oneexpd = 1.0-exp(-decay);

  unsigned int taildeg = 0, headdeg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *id=ML_IN_DEG(ll), *od=ML_OUT_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    taildeg += od[lt] + id[lt];
    headdeg += od[lh] + id[lh];
  }

  /* *** don't forget tail -> head */    
  CHANGE_STAT[0] =
    exp(decay) * (
		  ((1-pow(oneexpd,taildeg+degchange)) - (1-pow(oneexpd,taildeg))) +
		  ((1-pow(oneexpd,headdeg+degchange)) - (1-pow(oneexpd,headdeg)))
		  );
}

/*****************
 changestat: d_gwdegree_by_attr
*****************/
C_CHANGESTAT_FN(c_gwdegree_by_attr_ML_sum) { 
  /*The inputparams are assumed to be set up as follows:
    The first value is the decay parameter (as in Hunter et al, JASA 200?)
    The next sequence of values is the nodal attributes, coded as integers
         from 1 through N_CHANGE_STATS
  */
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double decay = *(inputs++);
  double *attrs = inputs-1;
  double oneexpd = 1.0-exp(-decay);

  unsigned int taildeg = 0, headdeg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *id=ML_IN_DEG(ll), *od=ML_OUT_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    taildeg += od[lt] + id[lt];
    headdeg += od[lh] + id[lh];
  }

    
  /* *** don't forget tail -> head */    
  CHANGE_STAT[(int)attrs[tail]] =
    exp(decay) * ((1-pow(oneexpd,taildeg+degchange)) - (1-pow(oneexpd,taildeg)));
  
  CHANGE_STAT[(int)attrs[head]] =
    exp(decay) * ((1-pow(oneexpd,headdeg+degchange)) - (1-pow(oneexpd,headdeg)));
}

/*****************
 changestat: d_gwb2degree_ML_sum
*****************/
C_CHANGESTAT_FN(c_gwb2degree_ML_sum) { 
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double decay = *(inputs++);
  double oneexpd = 1.0-exp(-decay);

  unsigned int b2deg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *id=ML_IN_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    b2deg += id[lh];
  }

    
  /* *** don't forget tail -> head */    
  CHANGE_STAT[0] =
        exp(decay) * ((1-pow(oneexpd,b2deg+degchange)) - (1-pow(oneexpd,b2deg)));
}

/*****************
 changestat: d_gwb2degree_by_attr_ML_sum
*****************/
C_CHANGESTAT_FN(c_gwb2degree_by_attr_ML_sum) { 
  /*The inputparams are assumed to be set up as follows:
    The first value is the decay parameter (as in Hunter et al, JASA 200?)
    The next sequence of values is the nodal attributes, coded as integers
         from 1 through N_CHANGE_STATS
  */
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double decay = *(inputs++);
  double *attrs = inputs-1;
  double oneexpd = 1.0-exp(-decay);

  unsigned int b2deg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *id=ML_IN_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    b2deg += id[lh];
  }

    
  /* *** don't forget tail -> head */    
  CHANGE_STAT[(int)attrs[head]] =
    exp(decay) * ((1-pow(oneexpd,b2deg+degchange)) - (1-pow(oneexpd,b2deg)));
}

/*****************
 changestat: d_gwidegree_ML_sum
*****************/
C_CHANGESTAT_FN(c_gwidegree_ML_sum) { 
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double decay = *(inputs++);
  double oneexpd = 1.0-exp(-decay);

  unsigned int headdeg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *id=ML_IN_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    headdeg += id[lh];
  }

    
  /* *** don't forget tail -> head */
  CHANGE_STAT[0] =
    exp(decay) * ((1-pow(oneexpd,headdeg+degchange)) - (1-pow(oneexpd,headdeg)));      
}

/*****************
 changestat: d_gwidegree_by_attr_ML_sum
*****************/
C_CHANGESTAT_FN(c_gwidegree_by_attr_ML_sum) { 
  /*The inputparams are assumed to be set up as follows:
    The first value is the decay parameter (as in Hunter et al, JASA 200?)
    The next sequence of values is the nodal attributes, coded as integers
         from 1 through N_CHANGE_STATS
  */
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double decay = *(inputs++);
  double *attrs = inputs-1;
  double oneexpd = 1.0-exp(-decay);

  unsigned int headdeg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *id=ML_IN_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    headdeg += id[lh];
  }

    
  /* *** don't forget tail -> head */    
  CHANGE_STAT[(int)attrs[head]] =
    exp(decay) * ((1-pow(oneexpd,headdeg+degchange)) - (1-pow(oneexpd,headdeg)));      
}

/*****************
 changestat: d_gwodegree_ML_sum
*****************/
C_CHANGESTAT_FN(c_gwodegree_ML_sum) { 
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double decay = *(inputs++);
  double oneexpd = 1.0-exp(-decay);

  unsigned int taildeg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *od=ML_OUT_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    taildeg += od[lt];
  }
    
  /* *** don't forget tail -> head */    
  CHANGE_STAT[0] =
    exp(decay) * ((1-pow(oneexpd,taildeg+degchange)) - (1-pow(oneexpd,taildeg)));      
}

/*****************
 changestat: d_gwodegree_by_attr_ML_sum
*****************/
C_CHANGESTAT_FN(c_gwodegree_by_attr_ML_sum) { 
  /*The inputparams are assumed to be set up as follows:
    The first value is the decay parameter (as in Hunter et al, JASA 200?)
    The next sequence of values is the nodal attributes, coded as integers
         from 1 through N_CHANGE_STATS
  */
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double decay = *(inputs++);
  double *attrs = inputs-1;
  double oneexpd = 1.0-exp(-decay);

  unsigned int taildeg = 0, headdeg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *od=ML_OUT_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    taildeg += od[lt];
  }

    
  /* *** don't forget tail -> head */    
  CHANGE_STAT[(int)attrs[tail]] =
    exp(decay) * ((1-pow(oneexpd,taildeg+degchange)) - (1-pow(oneexpd,taildeg)));
}

/*****************
 changestat: d_idegree_ML_sum
*****************/
C_CHANGESTAT_FN(c_idegree_ML_sum) { 
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double *degs = inputs;

  /* *** don't forget tail -> head */

  unsigned int headdeg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *id=ML_IN_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    headdeg += id[lh];
  }
  
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
    Vertex deg = (Vertex)degs[j];
    CHANGE_STAT[j] += (headdeg + degchange == deg) - (headdeg == deg);
  }
}

/*****************
 changestat: d_degree_by_attr_ML_sum
*****************/
C_CHANGESTAT_FN(c_idegree_by_attr_ML_sum) { 
  /* The inputparams are assumed to be set up as follows:
  The first 2*nstats values are in pairs:  (degree, attrvalue)
  The values following the first 2*nstats values are the nodal attributes.
  */
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);

  /* *** don't forget tail -> head */    

  unsigned int headdeg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *id=ML_IN_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    headdeg += id[lh];
  }   

  int headattr = inputs[2*N_CHANGE_STATS + head - 1]; 
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
    Vertex d = (Vertex)inputs[2*j];
    int testattr = inputs[2*j + 1]; 
    if (headattr == testattr)  /* we have head attr match */
      CHANGE_STAT[j] += (headdeg + degchange == d) - (headdeg == d);
  }
}

/*****************
 changestat: d_mutual_ML

 (1,1) -> anything = -1
 anything -> (1,1) = +1
*****************/
C_CHANGESTAT_FN(c_mutual_ML){
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll1, 0);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll2, 1);
  
  double matchval;
  int j, ninputs, noattr;

  ninputs = N_INPUT_PARAMS - N_NODES - 2;
  noattr = (N_INPUT_PARAMS == 2);

  /* *** don't forget tail -> head */
  Vertex lt = ML_IO_TAIL(ll1, tail), lh = ML_IO_HEAD(ll1, head);
  int l1th = ergm_LayerLogic2(lt, lh, tail, head, ll1, 2);
  int l1ht = ergm_LayerLogic2(lh, lt, tail, head, ll1, 2);
  int l2th = ergm_LayerLogic2(lt, lh, tail, head, ll2, 2);
  int l2ht = ergm_LayerLogic2(lh, lt, tail, head, ll2, 2);

  int change =
    +((l1th&2)&&(l2ht&2))-((l1th&1)&&(l2ht&1)) // t-l1->h and h->l2->t
    +((l2th&2)&&(l1ht&2))-((l2th&1)&&(l1ht&1)) // t-l2->h and h->l1->t
    ;
  
  if(change) { /* otherwise, no change occurs */
      if (noattr) { /* "plain vanilla" mutual, without node attributes */
        CHANGE_STAT[0] += change;
      } else { /* Only consider mutuals where node attributes match */
        matchval = INPUT_PARAM[tail+ninputs-1+2];
        if (matchval == INPUT_PARAM[head+ninputs-1+2]) { /* We have a match! */
          if (ninputs==0) {/* diff=F in network statistic specification */
            CHANGE_STAT[0] += change;
          } else { /* diff=T */
            for (j=0; j<ninputs; j++) {
              if (matchval == INPUT_PARAM[j+2]) 
                CHANGE_STAT[j] += change;
            }
          }
        }
      }
    }
}

/*****************
 changestat: d_mutual_by_attr_ML
*****************/
C_CHANGESTAT_FN(c_mutual_by_attr_ML) { 
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll1, 0);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll2, 1);

  int j, ninputs;

  ninputs = N_INPUT_PARAMS - N_NODES - 2;

  /* *** don't forget tail -> head */
  Vertex lt = ML_IO_TAIL(ll1, tail), lh = ML_IO_HEAD(ll1, head);
  int l1th = ergm_LayerLogic2(lt, lh, tail, head, ll1, 2);
  int l1ht = ergm_LayerLogic2(lh, lt, tail, head, ll1, 2);
  int l2th = ergm_LayerLogic2(lt, lh, tail, head, ll2, 2);
  int l2ht = ergm_LayerLogic2(lh, lt, tail, head, ll2, 2);

  int change =
    +((l1th&2)&&(l2ht&2))-((l1th&1)&&(l2ht&1)) // t-l1->h and h->l2->t
    +((l2th&2)&&(l1ht&2))-((l2th&1)&&(l1ht&1)) // t-l2->h and h->l1->t
    ;

    if (change) { /* otherwise, no change occurs */
      for (j=0; j<ninputs; j++) {
        if (INPUT_PARAM[tail+ninputs-1+2] == INPUT_PARAM[j+2]){CHANGE_STAT[j] += change;}
        if (INPUT_PARAM[head+ninputs-1+2] == INPUT_PARAM[j+2]){CHANGE_STAT[j] += change;}
      }
    }
}

/*****************
 changestat: d_odegree_ML_sum
*****************/
C_CHANGESTAT_FN(c_odegree_ML_sum) { 
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);
  double *degs = inputs;

  /* *** don't forget tail -> head */

  unsigned int taildeg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *od=ML_OUT_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    taildeg += od[lt];
  }
  
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
    Vertex deg = (Vertex)degs[j];
    CHANGE_STAT[j] += (taildeg + degchange == deg) - (taildeg == deg);
  }
}

/*****************
 changestat: d_odegree_by_attr_ML_sum
*****************/
C_CHANGESTAT_FN(c_odegree_by_attr_ML_sum) { 
  /* The inputparams are assumed to be set up as follows:
  The first 2*nstats values are in pairs:  (degree, attrvalue)
  The values following the first 2*nstats values are the nodal attributes.
  */
  double *inputs = INPUT_ATTRIB; // Already shifted past the auxiliaries.
  unsigned int nml = *(inputs++);

  /* *** don't forget tail -> head */    

  unsigned int taildeg = 0;
  int degchange = 0;

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    Vertex *od=ML_OUT_DEG(ll);
    Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head);
    degchange += ergm_LayerLogic(tail, head, ll, TRUE);
    taildeg += od[lt];
  }   

  int tailattr = inputs[2*N_CHANGE_STATS + tail - 1]; 
  for(unsigned int j = 0; j < N_CHANGE_STATS; j++) {
    Vertex d = (Vertex)inputs[2*j];
    int testattr = inputs[2*j + 1]; 
    if (tailattr == testattr)  /* we have tail attr match */
      CHANGE_STAT[j] += (taildeg + degchange == d) - (taildeg == d);
  }
}

