#include "changestat_operator.h"
#include "model.h"
#define STRICT_Wt_HEADERS
#include "wtchangestat.h"
#include "wtmodel.h"

#include "storage.h"

typedef struct{Network nw; Model *m;} StoreNetAndModel;

/* import_binary_term_sum 

   A term to wrap dyad-independent binary ergm terms by taking their
   change statistic from an empty network (i.e., their equivalent
   dyadic covariate value) and multiplying it by the difference
   between the previous and the new dyad value. 

*/

WtI_CHANGESTAT_FN(i_import_binary_term_sum){
  double *inputs = INPUT_PARAM;
  ALLOC_STORAGE(1, StoreNetAndModel, store);

  Model *m = store->m = unpack_Modelasdouble(&inputs);

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

  UPDATE_STORAGE(tail, head, m, mynwp);
}

WtF_CHANGESTAT_FN(f_import_binary_term_sum){
  GET_STORAGE(StoreNetAndModel, store);
  Model *m = store->m;
  Network *mynwp = &(store->nw);

  ModelDestroy(m, mynwp);
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

  mtp->storage = m = unpack_Modelasdouble(&inputs);

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
    UPDATE_STORAGE(tail, head, m, bnwp);
  }
}

WtF_CHANGESTAT_FN(f_import_binary_term_nonzero){
  GET_AUX_STORAGE(Network, bnwp);
  GET_STORAGE(Model, m);

  ModelDestroy(m, bnwp);
  mtp->storage = NULL;
}

/* _binary_nonzero_net 

   Maintain a binary network that mirrors the valued one in that it
   has an edge wherever the value is not 0.

*/

WtI_CHANGESTAT_FN(i__binary_nonzero_net){
  ALLOC_AUX_STORAGE(1, Network, bnwp);

  *bnwp = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, FALSE, 0, NULL);
  // FIXME: This is suboptimal, since all trees will be highly
  // unbalanced.
  WtEXEC_THROUGH_NET_EDGES(t, h, e, w, {
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
