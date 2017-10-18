#include "wtchangestats_operator.h"

/* import_binary_term_sum 

   A term to wrap dyad-independent binary ergm terms by taking their
   change statistic from an empty network (i.e., their equivalent
   dyadic covariate value) and multiplying it by the difference
   between the previous and the new dyad value. 

*/

WtI_CHANGESTAT_FN(i_import_binary_term_sum){
  double *inputs = INPUT_PARAM;
  ALLOC_STORAGE(1, StoreNetAndModel, store);

  Model *m = store->m = unpack_Model_as_double(&inputs);

  store->nw = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, FALSE, 0, NULL);
  Network *mynwp = &(store->nw);
  
  InitStats(mynwp, m);
}

WtC_CHANGESTAT_FN(c_import_binary_term_sum){
  GET_STORAGE(StoreNetAndModel, store);
  Model *m = store->m;
  Network *mynwp = &(store->nw);
  double oldweight = WtGETWT(tail,head);
    
  ChangeStats(1, &tail, &head, mynwp, m);

  for(unsigned int i=0; i<N_CHANGE_STATS; i++)
    CHANGE_STAT[i] = m->workspace[i]*(weight-oldweight);
}

WtU_CHANGESTAT_FN(u_import_binary_term_sum){
  GET_STORAGE(StoreNetAndModel, store);
  Model *m = store->m;
  Network *mynwp = &(store->nw);

  UPDATE_STORAGE(tail, head, mynwp, m, NULL);
}

WtF_CHANGESTAT_FN(f_import_binary_term_sum){
  GET_STORAGE(StoreNetAndModel, store);
  Model *m = store->m;
  Network *mynwp = &(store->nw);

  ModelDestroy(mynwp, m);
  NetworkDestroy(mynwp);
}

/* import_binary_term_nonzero 

   A term to wrap abitrary binary ergm terms by constructing a binary
   network that mirrors the valued one in that it has an edge wherever
   the value is not 0.

*/

WtI_CHANGESTAT_FN(i_import_binary_term_nonzero){
  double *inputs = INPUT_PARAM;
  GET_AUX_STORAGE(Network, bnwp); inputs++;
  GET_STORAGE(Model, m); // Only need the pointer, no allocation needed.

  STORAGE = m = unpack_Model_as_double(&inputs);

  InitStats(bnwp, m);
}

WtC_CHANGESTAT_FN(c_import_binary_term_nonzero){
  GET_AUX_STORAGE(Network, bnwp);
  GET_STORAGE(Model, m);
  
  double oldweight = WtGETWT(tail,head);

  if((weight!=0)!=(oldweight!=0)){ // If going from 0 to nonzero or vice versa...
    ChangeStats(1, &tail, &head, bnwp, m);
  }
  
  memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double));
}

WtU_CHANGESTAT_FN(u_import_binary_term_nonzero){
  GET_AUX_STORAGE(Network, bnwp);
  GET_STORAGE(Model, m);
  
  double oldweight = WtGETWT(tail,head);

  if((weight!=0)!=(oldweight!=0)){ // If going from 0 to nonzero or vice versa...
    UPDATE_STORAGE(tail, head, bnwp, m, NULL);
  }
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
  double *inputs = INPUT_PARAM;
  GET_AUX_STORAGE(StoreNetAndWtModel, storage); inputs++;
  Network *bnwp = &(storage->nw);

  GET_STORAGE(Model, m); // Only need the pointer, no allocation needed.

  STORAGE = m = unpack_Model_as_double(&inputs);

  InitStats(bnwp, m);
}

WtC_CHANGESTAT_FN(c_import_binary_term_form){
  GET_AUX_STORAGE(StoreNetAndWtModel, storage);
  Network *bnwp = &(storage->nw);
  GET_STORAGE(Model, m);

  WtChangeStats(1, &tail, &head, &weight, nwp, storage->m);
  
  if(*(storage->m->workspace)!=0){ // If the binary view changes...
    ChangeStats(1, &tail, &head, bnwp, m);
    memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double));
  } // Otherwise, leave the change stats at 0.
}

WtU_CHANGESTAT_FN(u_import_binary_term_form){
  GET_AUX_STORAGE(StoreNetAndWtModel, storage);
  Network *bnwp = &(storage->nw);
  GET_STORAGE(Model, m);

  WtChangeStats(1, &tail, &head, &weight, nwp, storage->m);
  
  if(*(storage->m->workspace)!=0){ // If the binary view changes...
    UPDATE_STORAGE(tail, head, bnwp, m, NULL);
  }
}

WtF_CHANGESTAT_FN(f_import_binary_term_form){
  GET_AUX_STORAGE(StoreNetAndWtModel, storage);
  Network *bnwp = &(storage->nw);
  GET_STORAGE(Model, m);

  ModelDestroy(bnwp, m);
  STORAGE = NULL;
}

/* _binary_nonzero_net 

   Maintain a binary network that mirrors the valued one in that it
   has an edge wherever the value is not 0.

*/

WtI_CHANGESTAT_FN(i__binary_nonzero_net){
  ALLOC_AUX_STORAGE(1, Network, bnwp);

  *bnwp = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, FALSE, 0, NULL);
  WtEXEC_THROUGH_NET_EDGES_PRE(t, h, e, w, {
      if(w!=0) ToggleEdge(t, h, bnwp);
    });
}

WtU_CHANGESTAT_FN(u__binary_nonzero_net){
  GET_AUX_STORAGE(Network, bnwp);
  double oldweight = WtGETWT(tail,head);

  if((weight!=0)!=(oldweight!=0)){ // If going from 0 to nonzero or vice versa...
    ToggleEdge(tail, head, bnwp);
  }
}

WtF_CHANGESTAT_FN(f__binary_nonzero_net){
  GET_AUX_STORAGE(Network, bnwp);
  NetworkDestroy(bnwp);
}


/* _binary_formula_net 

   Maintain a binary network that mirrors the valued one in that it
   has an edge wherever the contribution of the a given term (edges,
   nonzero, ininterval, atleast, atmost, etc.) whose dyadwise value is
   either 0 or 1 is 1.

*/


WtI_CHANGESTAT_FN(i__binary_formula_net){
  double *inputs = INPUT_PARAM;
  ALLOC_AUX_STORAGE(1, StoreNetAndWtModel, storage); inputs++;
  WtModel *m = storage->m = unpack_WtModel_as_double(&inputs);
  Network *bnwp = &(storage->nw);
  
  WtInitStats(nwp, m);
  
  *bnwp = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, FALSE, 0, NULL);
  WtEXEC_THROUGH_NET_EDGES_PRE(t, h, e, w, {
      if(w!=0){
	WtChangeStats(1, &t, &h, 0, nwp, m);
	// I.e., if reducing the value from the current value to 0
	// decreases the statistic, add edge to the binary network.
	if(*(m->workspace)==-1) 
	  AddEdgeToTrees(t, h, bnwp);
	else if(*(m->workspace)!=0) error("Binary test term may have a dyadwise contribution of either 0 or 1. Memory has not been deallocated, so restart R soon.");
      }
    });
}

WtU_CHANGESTAT_FN(u__binary_formula_net){
  GET_AUX_STORAGE(StoreNetAndWtModel, storage);
  WtModel *m = storage->m;
  Network *bnwp = &(storage->nw);

  WtChangeStats(1, &tail, &head, &weight, nwp, m);
  switch((int) *(m->workspace)){
  case  0: break;
  case -1: DeleteEdgeFromTrees(tail,head,bnwp); break;
  case +1: AddEdgeToTrees(tail,head,bnwp); break;
  default: error("Binary test term may have a dyadwise contribution of either 0 or 1. Memory has not been deallocated, so restart R soon."); 
  }

  WtUPDATE_STORAGE(tail, head, weight, nwp, m, NULL);
}

WtF_CHANGESTAT_FN(f__binary_formula_net){
  GET_AUX_STORAGE(StoreNetAndWtModel, storage);
  WtModel *m = storage->m;
  Network *bnwp = &(storage->nw);
  WtModelDestroy(nwp, m);
  NetworkDestroy(bnwp);
  // WtDestroyStats() will deallocate the rest.
}
