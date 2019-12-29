#include "wtchangestats_operator.h"

/* passthrough(formula) */

WtI_CHANGESTAT_FN(i_wtpassthrough_term){
  // No need to allocate it: we are only storing a pointer to a model.
  WtModel *m = STORAGE = WtModelInitialize(getListElement(mtp->R, "submodel"), NULL,  nwp, FALSE);

  WtSELECT_C_OR_D_BASED_ON_SUBMODEL(m);
  WtDELETE_IF_UNUSED_IN_SUBMODEL(u_func, m);
  WtDELETE_IF_UNUSED_IN_SUBMODEL(z_func, m);
}

WtD_CHANGESTAT_FN(d_wtpassthrough_term){
  GET_STORAGE(WtModel, m);

  WtChangeStats(ntoggles, tails, heads, weights, nwp, m);

  memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double));
}

WtC_CHANGESTAT_FN(c_wtpassthrough_term){
  GET_STORAGE(WtModel, m);

  WtChangeStats1(tail, head, weight, nwp, m, edgeweight);

  memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double));
}

WtU_CHANGESTAT_FN(u_wtpassthrough_term){
  GET_STORAGE(WtModel, m);

  WtUPDATE_STORAGE(tail, head, weight, nwp, m, NULL, edgeweight);
}

WtZ_CHANGESTAT_FN(z_wtpassthrough_term){
  GET_STORAGE(WtModel, m);

  WtZStats(nwp, m);

  memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double));
}

WtF_CHANGESTAT_FN(f_wtpassthrough_term){
  GET_STORAGE(WtModel, m);

  WtModelDestroy(nwp, m);

  STORAGE = NULL;
}

/* import_binary_term_sum 

   A term to wrap dyad-independent binary ergm terms by taking their
   change statistic from an empty network (i.e., their equivalent
   dyadic covariate value) and multiplying it by the difference
   between the previous and the new dyad value. 

*/

WtI_CHANGESTAT_FN(i_import_binary_term_sum){
  ALLOC_STORAGE(1, StoreNetAndModel, store);

  store->nwp = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, FALSE, 0, NULL);
  Network *mynwp = store->nwp;
  store->m = ModelInitialize(getListElement(mtp->R, "submodel"), NULL,  mynwp, FALSE);
  DELETE_IF_UNUSED_IN_SUBMODEL(u_func, store->m);
  DELETE_IF_UNUSED_IN_SUBMODEL(z_func, store->m);
}

WtC_CHANGESTAT_FN(c_import_binary_term_sum){
  GET_STORAGE(StoreNetAndModel, store);
  Model *m = store->m;
  Network *mynwp = store->nwp;
    
  ChangeStats1(tail, head, mynwp, m, FALSE); // mynwp is a dummy network that is always empty.

  for(unsigned int i=0; i<N_CHANGE_STATS; i++)
    CHANGE_STAT[i] = m->workspace[i]*(weight-edgeweight);
}

WtU_CHANGESTAT_FN(u_import_binary_term_sum){
  GET_STORAGE(StoreNetAndModel, store);
  Model *m = store->m;
  Network *mynwp = store->nwp;

  GET_EDGE_UPDATE_STORAGE(tail, head, mynwp, m, NULL);
}

/* WtZ_CHANGESTAT_FN(z_import_binary_term_sum) is not meaningful. */

WtF_CHANGESTAT_FN(f_import_binary_term_sum){
  GET_STORAGE(StoreNetAndModel, store);
  Model *m = store->m;
  Network *mynwp = store->nwp;

  ModelDestroy(mynwp, m);
  NetworkDestroy(mynwp);
}

/* import_binary_term_nonzero 

   A term to wrap abitrary binary ergm terms by constructing a binary
   network that mirrors the valued one in that it has an edge wherever
   the value is not 0.

*/

WtI_CHANGESTAT_FN(i_import_binary_term_nonzero){
  GET_AUX_STORAGE(Network, bnwp);
  GET_STORAGE(Model, m); // Only need the pointer, no allocation needed.

  STORAGE = m = ModelInitialize(getListElement(mtp->R, "submodel"), NULL,  bnwp, FALSE);
  DELETE_IF_UNUSED_IN_SUBMODEL(u_func, m);
  DELETE_IF_UNUSED_IN_SUBMODEL(z_func, m);
}

WtC_CHANGESTAT_FN(c_import_binary_term_nonzero){
  GET_AUX_STORAGE(Network, bnwp);
  GET_STORAGE(Model, m);
  

  if((weight!=0)!=(edgeweight!=0)){ // If going from 0 to nonzero or vice versa...
    ChangeStats1(tail, head, bnwp, m, edgeweight!=0);
  }
  
  memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double));
}

WtU_CHANGESTAT_FN(u_import_binary_term_nonzero){
  GET_AUX_STORAGE(Network, bnwp);
  GET_STORAGE(Model, m);

  if((weight!=0)!=(edgeweight!=0)){ // If going from 0 to nonzero or vice versa...
    GET_EDGE_UPDATE_STORAGE(tail, head, bnwp, m, NULL);
  }
}

WtZ_CHANGESTAT_FN(z_import_binary_term_nonzero){
  GET_AUX_STORAGE(Network, bnwp);
  GET_STORAGE(Model, m);

  ZStats(bnwp, m);

  memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double));
}

WtF_CHANGESTAT_FN(f_import_binary_term_nonzero){
  GET_AUX_STORAGE(Network, bnwp);
  GET_STORAGE(Model, m);

  ModelDestroy(bnwp, m);
  STORAGE = NULL;
}


/* import_binary_term_form 

   A term to wrap abitrary binary ergm terms by constructing a binary
   network that mirrors the valued one in that it has an edge iff the term
   in the second formula contributes +1 due to that dyad.

*/

WtI_CHANGESTAT_FN(i_import_binary_term_form){
  GET_AUX_STORAGE(StoreNetAndWtModel, storage);
  Network *bnwp = storage->nwp;

  GET_STORAGE(Model, m); // Only need the pointer, no allocation needed.

  STORAGE = m = ModelInitialize(getListElement(mtp->R, "submodel"), NULL, bnwp, FALSE);
  DELETE_IF_UNUSED_IN_SUBMODEL(u_func, m);
  DELETE_IF_UNUSED_IN_SUBMODEL(z_func, m);
}

WtC_CHANGESTAT_FN(c_import_binary_term_form){
  GET_AUX_STORAGE(StoreNetAndWtModel, storage);
  Network *bnwp = storage->nwp;
  GET_STORAGE(Model, m);

  WtChangeStats1(tail, head, weight, nwp, storage->m, edgeweight);
  
  if(*(storage->m->workspace)!=0){ // If the binary view changes...
    ChangeStats1(tail, head, bnwp, m, IS_OUTEDGE(tail, head, bnwp));
    memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double));
  } // Otherwise, leave the change stats at 0.
}

WtU_CHANGESTAT_FN(u_import_binary_term_form){
  GET_AUX_STORAGE(StoreNetAndWtModel, storage);
  Network *bnwp = storage->nwp;
  GET_STORAGE(Model, m);

  WtChangeStats1(tail, head, weight, nwp, storage->m, edgeweight);
  
  if(*(storage->m->workspace)!=0){ // If the binary view changes...
    GET_EDGE_UPDATE_STORAGE(tail, head, bnwp, m, NULL);
  }
}

WtZ_CHANGESTAT_FN(z_import_binary_term_form){
  GET_AUX_STORAGE(StoreNetAndWtModel, storage);
  Network *bnwp = storage->nwp;
  GET_STORAGE(Model, m);

  ZStats(bnwp, m);
  memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double));
}

WtF_CHANGESTAT_FN(f_import_binary_term_form){
  GET_AUX_STORAGE(StoreNetAndWtModel, storage);
  Network *bnwp = storage->nwp;
  GET_STORAGE(Model, m);

  ModelDestroy(bnwp, m);
  STORAGE = NULL;
}

/* _binary_nonzero_net 

   Maintain a binary network that mirrors the valued one in that it
   has an edge wherever the value is not 0.

*/

WtI_CHANGESTAT_FN(i__binary_nonzero_net){
  Network *bnwp = AUX_STORAGE = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, FALSE, 0, NULL);
  WtEXEC_THROUGH_NET_EDGES_PRE(t, h, e, w, {
      if(w!=0) ToggleEdge(t, h, bnwp);
    });
}

WtU_CHANGESTAT_FN(u__binary_nonzero_net){
  GET_AUX_STORAGE(Network, bnwp);

  if((weight!=0)!=(edgeweight!=0)){ // If going from 0 to nonzero or vice versa...
    ToggleEdge(tail, head, bnwp);
  }
}

WtF_CHANGESTAT_FN(f__binary_nonzero_net){
  GET_AUX_STORAGE(Network, bnwp);
  NetworkDestroy(bnwp);
  AUX_STORAGE = NULL;
}


/* _binary_formula_net 

   Maintain a binary network that mirrors the valued one in that it
   has an edge wherever the contribution of the a given term (edges,
   nonzero, ininterval, atleast, atmost, etc.) whose dyadwise value is
   either 0 or 1 is 1.

*/


WtI_CHANGESTAT_FN(i__binary_formula_net){
  ALLOC_AUX_STORAGE(1, StoreNetAndWtModel, storage);
  WtModel *m = storage->m = WtModelInitialize(getListElement(mtp->R, "submodel"), NULL, nwp, FALSE);
  Network *bnwp = storage->nwp = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, FALSE, 0, NULL);
 
  WtEXEC_THROUGH_NET_EDGES_PRE(t, h, e, w, {
      if(w!=0){
	WtChangeStats1(t, h, 0, nwp, m, w);
	// I.e., if reducing the value from the current value to 0
	// decreases the statistic, add edge to the binary network.
	if(*(m->workspace)==-1) 
	  AddEdgeToTrees(t, h, bnwp);
	else if(*(m->workspace)!=0) error("Binary test term may have a dyadwise contribution of either 0 or 1. Memory has not been deallocated, so restart R soon.");
      }
    });

  WtDELETE_IF_UNUSED_IN_SUBMODEL(z_func, m);
}

WtU_CHANGESTAT_FN(u__binary_formula_net){
  GET_AUX_STORAGE(StoreNetAndWtModel, storage);
  WtModel *m = storage->m;
  Network *bnwp = storage->nwp;

  WtChangeStats1(tail, head, weight, nwp, m, edgeweight);
  switch((int) *(m->workspace)){
  case  0: break;
  case -1: DeleteEdgeFromTrees(tail,head,bnwp); break;
  case +1: AddEdgeToTrees(tail,head,bnwp); break;
  default: error("Binary test term may have a dyadwise contribution of either 0 or 1. Memory has not been deallocated, so restart R soon."); 
  }

  WtUPDATE_STORAGE(tail, head, weight, nwp, m, NULL, edgeweight);
}

WtF_CHANGESTAT_FN(f__binary_formula_net){
  GET_AUX_STORAGE(StoreNetAndWtModel, storage);
  WtModel *m = storage->m;
  Network *bnwp = storage->nwp;
  WtModelDestroy(nwp, m);
  NetworkDestroy(bnwp);
  // WtDestroyStats() will deallocate the rest.
}


// wtSum: Take a weighted sum of the models' statistics.

WtI_CHANGESTAT_FN(i_wtSum){
  /*
    inputs expected:
    1: number of models (nms)
    1: total length of all weight matrices (tml)
    tml: a list of mapping matrices in row-major order
    nms*?: submodel specifications for nms submodels
  */
  
  double *inputs = INPUT_PARAM; 
  unsigned int nms = *(inputs++);
  unsigned int tml = *(inputs++);
  inputs+=tml; // Jump to start of model specifications.
  
  ALLOC_STORAGE(nms, WtModel*, ms);

  SEXP submodels = getListElement(mtp->R, "submodels");
  for(unsigned int i=0; i<nms; i++){
    ms[i] = WtModelInitialize(VECTOR_ELT(submodels,i), NULL, nwp, FALSE);
  }
  WtDELETE_IF_UNUSED_IN_SUBMODELS(u_func, ms, nms);
  WtDELETE_IF_UNUSED_IN_SUBMODELS(z_func, ms, nms);
}

WtC_CHANGESTAT_FN(c_wtSum){
  double *inputs = INPUT_PARAM; 
  GET_STORAGE(WtModel*, ms);
  unsigned int nms = *(inputs++);
  inputs++; //  Skip total length of weight matrices.
  double *wts = inputs;

  for(unsigned int i=0; i<nms; i++){
    WtModel *m = ms[i];
    WtChangeStats1(tail, head, weight, nwp, m, edgeweight);
    for(unsigned int j=0; j<m->n_stats; j++)
      for(unsigned int k=0; k<N_CHANGE_STATS; k++)
	CHANGE_STAT[k] += m->workspace[j]* *(wts++);
  }
}

WtU_CHANGESTAT_FN(u_wtSum){
  double *inputs = INPUT_PARAM; 
  GET_STORAGE(WtModel*, ms);
  unsigned int nms = *(inputs++);

  for(unsigned int i=0; i<nms; i++){
    WtModel *m = ms[i];
    WtUPDATE_STORAGE(tail, head, weight, nwp, m, NULL, edgeweight);
  }
}

WtZ_CHANGESTAT_FN(z_wtSum){
  double *inputs = INPUT_PARAM;
  GET_STORAGE(WtModel*, ms);
  unsigned int nms = *(inputs++);
  inputs++; //  Skip total length of weight matrices.
  double *wts = inputs;

  for(unsigned int i=0; i<nms; i++){
    WtModel *m = ms[i];
    WtZStats(nwp, m);
    for(unsigned int j=0; j<m->n_stats; j++)
      for(unsigned int k=0; k<N_CHANGE_STATS; k++)
	CHANGE_STAT[k] += m->workspace[j]* *(wts++);
  }
}

WtF_CHANGESTAT_FN(f_wtSum){
  double *inputs = INPUT_PARAM; 
  GET_STORAGE(WtModel*, ms);
  unsigned int nms = *(inputs++);

  for(unsigned int i=0; i<nms; i++){
    WtModelDestroy(nwp, ms[i]);
  }
}
