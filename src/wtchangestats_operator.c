
#include "changestat.h"
#include "model.h"
#include "wtchangestat.h"
#include "wtmodel.h"
#include "storage.h"

typedef struct{Network nw; Model *m;} StoreNetAndModel;

static inline unsigned char *unpack_doubletochar(double **x){
  unsigned int l = (*x)[0];
  unsigned char *s = (unsigned char *) malloc((l+1)*sizeof(char));
  for(unsigned int i=0; i<l; i++){
    s[i] = (unsigned char) (unsigned int) (*x)[i+1];
  }
  s[l] = (unsigned char) 0;
  (*x)+=l+1;
  return s;
}

/* import_binary_term_sum 

   A term to wrap dyad-independent binary ergm terms by taking their
   change statistic from an empty network (i.e., their equivalent
   dyadic covariate value) and multiplying it by the difference
   between the previous and the new dyad value. 

*/

WtI_CHANGESTAT_FN(i_import_binary_term_sum){
  double *inputs = INPUT_PARAM;
  ALLOC_STORAGE(1, StoreNetAndModel, store);

  int n_terms = *(inputs++);
  char *fnames = (char *) unpack_doubletochar(&inputs);
  char *snames = (char *) unpack_doubletochar(&inputs);
  Model *m = store->m = ModelInitialize(fnames, snames, &inputs, n_terms);
  free(fnames);
  free(snames);

  store->nw = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, FALSE, 0, NULL);
  Network *mynwp = &(store->nw);
  
  InitStats(mynwp, m);
}

WtC_CHANGESTAT_FN(c_import_binary_term_sum){
  GET_STORAGE(StoreNetAndModel, store);
  Model *m = store->m;
  Network *mynwp = &(store->nw);
  double oldweight = GETWT(tail,head);
    
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

   FIXME: Use auxiliary storage so that multiple terms can share the network.

*/


WtI_CHANGESTAT_FN(i_import_binary_term_nonzero){
  double *inputs = INPUT_PARAM;
  ALLOC_STORAGE(1, StoreNetAndModel, store);

  int n_terms = *(inputs++);
  char *fnames = (char *) unpack_doubletochar(&inputs);
  char *snames = (char *) unpack_doubletochar(&inputs);
  Model *m = store->m = ModelInitialize(fnames, snames, &inputs, n_terms);
  free(fnames);
  free(snames);

  store->nw = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, FALSE, 0, NULL);
  Network *mynwp = &(store->nw);
  // FIXME: This is suboptimal, since all trees will be highly
  // unbalanced.
  EXEC_THROUGH_NET_EDGES(t, h, e, w, {
      if(w!=0) ToggleEdge(t, h, mynwp);
    });
  
  InitStats(mynwp, m);
}

WtC_CHANGESTAT_FN(c_import_binary_term_nonzero){
  GET_STORAGE(StoreNetAndModel, store);
  Model *m = store->m;
  Network *mynwp = &(store->nw);
  double oldweight = GETWT(tail,head);

  if((weight!=0)!=(oldweight!=0)){ // If going from 0 to nonzero or vice versa...
    ChangeStats(1, &tail, &head, mynwp, m);

    for(unsigned int i=0; i<N_CHANGE_STATS; i++)
      CHANGE_STAT[i] = m->workspace[i];
  }
}

WtU_CHANGESTAT_FN(u_import_binary_term_nonzero){
  GET_STORAGE(StoreNetAndModel, store);
  Model *m = store->m;
  Network *mynwp = &(store->nw);
  double oldweight = GETWT(tail,head);

  if((weight!=0)!=(oldweight!=0)){ // If going from 0 to nonzero or vice versa...
    UPDATE_STORAGE(tail, head, m, mynwp);
    ToggleEdge(tail, head, mynwp);
  }
}

WtF_CHANGESTAT_FN(f_import_binary_term_nonzero){
  GET_STORAGE(StoreNetAndModel, store);
  Model *m = store->m;
  Network *mynwp = &(store->nw);

  ModelDestroy(m, mynwp);
  NetworkDestroy(mynwp);
}
