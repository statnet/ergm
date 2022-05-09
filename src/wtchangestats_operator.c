/*  File src/wtchangestats_operator.c in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2022 Statnet Commons
 */
#include "wtchangestats_operator.h"
#include "ergm_wtchangestats_operator.h"
#include "ergm_util.h"

/* passthrough(formula) */

WtI_CHANGESTAT_FN(i_wtpassthrough_term){
  // No need to allocate it: we are only storing a pointer to a model.
  WtModel *m = STORAGE = WtModelInitialize(getListElement(mtp->R, "submodel"), mtp->ext_state,  nwp, FALSE);

  WtSELECT_C_OR_D_BASED_ON_SUBMODEL(m);
  WtDELETE_IF_UNUSED_IN_SUBMODEL(x_func, m);
  WtDELETE_IF_UNUSED_IN_SUBMODEL(z_func, m);
}

WtD_CHANGESTAT_FN(d_wtpassthrough_term){
  GET_STORAGE(WtModel, m);

  WtChangeStats(ntoggles, tails, heads, weights, nwp, m);

  memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double));
}

WtC_CHANGESTAT_FN(c_wtpassthrough_term){
  GET_STORAGE(WtModel, m);

  WtChangeStats1(tail, head, weight, nwp, m, edgestate);

  memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double));
}

WtZ_CHANGESTAT_FN(z_wtpassthrough_term){
  GET_STORAGE(WtModel, m);

  WtZStats(nwp, m, FALSE);

  memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double));
}

WtX_CHANGESTAT_PROPAGATE_FN(x_wtpassthrough_term, GET_STORAGE(WtModel, m), m)

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
  DELETE_IF_UNUSED_IN_SUBMODEL(z_func, store->m);
}

WtC_CHANGESTAT_FN(c_import_binary_term_sum){
  GET_STORAGE(StoreNetAndModel, store);
  Model *m = store->m;
  Network *mynwp = store->nwp;
    
  ChangeStats1(tail, head, mynwp, m, FALSE); // mynwp is a dummy network that is always empty.

  for(unsigned int i=0; i<N_CHANGE_STATS; i++)
    CHANGE_STAT[i] = m->workspace[i]*(weight-edgestate);
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
  DELETE_IF_UNUSED_IN_SUBMODEL(z_func, m);
}

WtC_CHANGESTAT_FN(c_import_binary_term_nonzero){
  GET_AUX_STORAGE(Network, bnwp);
  GET_STORAGE(Model, m);
  

  if((weight!=0)!=(edgestate!=0)){ // If going from 0 to nonzero or vice versa...
    ChangeStats1(tail, head, bnwp, m, edgestate!=0);
  }
  
  memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double));
}

WtZ_CHANGESTAT_FN(z_import_binary_term_nonzero){
  GET_AUX_STORAGE(Network, bnwp);
  GET_STORAGE(Model, m);

  ZStats(bnwp, m, FALSE);

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
  DELETE_IF_UNUSED_IN_SUBMODEL(z_func, m);
}

WtC_CHANGESTAT_FN(c_import_binary_term_form){
  GET_AUX_STORAGE(StoreNetAndWtModel, storage);
  Network *bnwp = storage->nwp;
  GET_STORAGE(Model, m);

  WtChangeStats1(tail, head, weight, nwp, storage->m, edgestate);
  
  if(*(storage->m->workspace)!=0){ // If the binary view changes...
    ChangeStats1(tail, head, bnwp, m, IS_OUTEDGE(tail, head, bnwp));
    memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double));
  } // Otherwise, leave the change stats at 0.
}

WtZ_CHANGESTAT_FN(z_import_binary_term_form){
  GET_AUX_STORAGE(StoreNetAndWtModel, storage);
  Network *bnwp = storage->nwp;
  GET_STORAGE(Model, m);

  ZStats(bnwp, m, FALSE);
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

  if((weight!=0)!=(edgestate!=0)){ // If going from 0 to nonzero or vice versa...
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

  WtChangeStats1(tail, head, weight, nwp, m, edgestate);
  switch((int) *(m->workspace)){
  case  0: break;
  case -1: DeleteEdgeFromTrees(tail,head,bnwp); break;
  case +1: AddEdgeToTrees(tail,head,bnwp); break;
  default: error("Binary test term may have a dyadwise contribution of either 0 or 1. Memory has not been deallocated, so restart R soon."); 
  }
}

WtF_CHANGESTAT_FN(f__binary_formula_net){
  GET_AUX_STORAGE(StoreNetAndWtModel, storage);
  WtModel *m = storage->m;
  Network *bnwp = storage->nwp;
  WtModelDestroy(nwp, m);
  NetworkDestroy(bnwp);
  // WtDestroyStats() will deallocate the rest.
}


/* .submodel_and_summary(formula) */

WtI_CHANGESTAT_FN(i__wtsubmodel_and_summary_term){
  ALLOC_AUX_STORAGE(1, StoreWtModelAndStats, storage);

  // Unpack the submodel.
  WtModel *m = storage->m = WtModelInitialize(getListElement(mtp->R, "submodel"),  NULL, nwp, FALSE);

  storage->stats = Calloc(m->n_stats, double);

  WtSummStats(0, NULL, NULL, NULL, nwp, m);
  memcpy(storage->stats, m->workspace, m->n_stats*sizeof(double));

  DELETE_IF_UNUSED_IN_SUBMODEL(z_func, m);
}

WtU_CHANGESTAT_FN(u__wtsubmodel_and_summary_term){
  GET_AUX_STORAGE(StoreWtModelAndStats, storage);
  WtModel *m = storage->m;

  WtChangeStats1(tail, head, weight, nwp, m, edgestate);
  addonto(storage->stats, m->workspace, m->n_stats);
}

WtF_CHANGESTAT_FN(f__wtsubmodel_and_summary_term){
  GET_AUX_STORAGE(StoreWtModelAndStats, storage);

  Free(storage->stats);
  WtModelDestroy(nwp, storage->m);
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
    ms[i] = WtModelInitialize(VECTOR_ELT(submodels,i), isNULL(mtp->ext_state) ? NULL : VECTOR_ELT(mtp->ext_state,i), nwp, FALSE);
  }
  WtDELETE_IF_UNUSED_IN_SUBMODELS(x_func, ms, nms);
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
    WtChangeStats1(tail, head, weight, nwp, m, edgestate);
    for(unsigned int j=0; j<m->n_stats; j++)
      for(unsigned int k=0; k<N_CHANGE_STATS; k++)
	CHANGE_STAT[k] += m->workspace[j]* *(wts++);
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
    WtZStats(nwp, m, FALSE);
    for(unsigned int j=0; j<m->n_stats; j++)
      for(unsigned int k=0; k<N_CHANGE_STATS; k++)
	CHANGE_STAT[k] += m->workspace[j]* *(wts++);
  }
}

WtX_CHANGESTAT_FN(x_wtSum){
  double *inputs = INPUT_PARAM;
  GET_STORAGE(WtModel*, ms);
  unsigned int nms = *(inputs++);
  inputs++; //  Skip total length of weight matrices.
  double *wts = inputs;

  for(unsigned int i=0; i<nms; i++){
    WtModel *m = ms[i];
    WtPROPAGATE_X_SIGNAL_INTO(nwp, m, m->workspace);
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


// Log: Take a natural logarithm of the model's statistics.

WtI_CHANGESTAT_FN(i_wtLog){
  GET_AUX_STORAGE(StoreWtModelAndStats, modstats);
  DELETE_IF_UNUSED_IN_SUBMODEL(z_func, modstats->m);
}

WtC_CHANGESTAT_FN(c_wtLog){
  GET_AUX_STORAGE(StoreWtModelAndStats, modstats);
  double *log0 = INPUT_PARAM;

  WtChangeStats1(tail, head, weight, nwp, modstats->m, edgestate);
  for(unsigned int i=0; i<N_CHANGE_STATS; i++){
    if(modstats->m->workspace[i] == 0) CHANGE_STAT[i] = 0;
    else{
      double old = modstats->stats[i];
      old = old == 0 ? log0[i] : log(old);
      double new = modstats->stats[i]+modstats->m->workspace[i];
      new = new == 0 ? log0[i] : log(new);
      CHANGE_STAT[i] = new - old;
    }
  }
}

WtZ_CHANGESTAT_FN(z_wtLog){
  GET_AUX_STORAGE(StoreWtModelAndStats, modstats);
  double* log0 = INPUT_PARAM;

  WtEmptyNetworkStats(modstats->m, FALSE);
  memcpy(CHANGE_STAT, modstats->m->workspace, N_CHANGE_STATS*sizeof(double));
  WtZStats(nwp, modstats->m, FALSE);
  for(unsigned int i=0; i<N_CHANGE_STATS; i++){
    if(modstats->m->workspace[i] == 0) CHANGE_STAT[i] = 0;
    else{
      double old = CHANGE_STAT[i];
      old = old == 0 ? log0[i] : log(old);
      double new = CHANGE_STAT[i]+modstats->m->workspace[i];
      new = new == 0 ? log0[i] : log(new);
      CHANGE_STAT[i] = new - old;
    }
  }
}

// Exp: Exponentiate the model's statistics.

WtI_CHANGESTAT_FN(i_wtExp){
  GET_AUX_STORAGE(StoreWtModelAndStats, modstats);
  DELETE_IF_UNUSED_IN_SUBMODEL(z_func, modstats->m);
}

WtC_CHANGESTAT_FN(c_wtExp){
  GET_AUX_STORAGE(StoreWtModelAndStats, modstats);

  WtChangeStats1(tail, head, weight, nwp, modstats->m, edgestate);
  for(unsigned int i=0; i<N_CHANGE_STATS; i++){
    if(modstats->m->workspace[i] == 0) CHANGE_STAT[i] = 0;
    else CHANGE_STAT[i] = exp(modstats->stats[i]+modstats->m->workspace[i]) - exp(modstats->stats[i]);
  }
}

WtZ_CHANGESTAT_FN(z_wtExp){
  GET_AUX_STORAGE(StoreWtModelAndStats, modstats);
  WtEmptyNetworkStats(modstats->m, FALSE);
  memcpy(CHANGE_STAT, modstats->m->workspace, N_CHANGE_STATS*sizeof(double));
  WtZStats(nwp, modstats->m, FALSE);
  for(unsigned int i=0; i<N_CHANGE_STATS; i++){
    if(modstats->m->workspace[i] == 0) CHANGE_STAT[i] = 0;
    else CHANGE_STAT[i] = exp(CHANGE_STAT[i]+modstats->m->workspace[i]) - exp(CHANGE_STAT[i]);
  }
}
