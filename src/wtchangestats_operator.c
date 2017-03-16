#include "wtchangestat.h"
#include "wtmodel.h"

#include "changestat.h"
#include "model.h"

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

WtI_CHANGESTAT_FN(i_import_binary_term_sum){
  double *inputs = INPUT_PARAM;
  ALLOC_STORAGE(1, StoreNetAndModel, store);

  int n_terms = *(inputs++);
  char *fnames = (char *) unpack_doubletochar(&inputs);
  char *snames = (char *) unpack_doubletochar(&inputs);
  Model *m = store->m = ModelInitialize(fnames, snames, &inputs, n_terms);
  free(fnames);
  free(snames);

  store->nw = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, FALSE, 0, NULL, m->n_aux);
  Network *mynwp = &(store->nw);
  
  InitStats(mynwp, m);
}

WtC_CHANGESTAT_FN(c_import_binary_term_sum){
  GET_STORAGE(StoreNetAndModel, store);
  Model *m = store->m;
  Network *mynwp = &(store->nw);
    
  ChangeStats(1, &tail, &head, mynwp, m);

  for(unsigned int i=0; i<N_CHANGE_STATS; i++)
    CHANGE_STAT[i] = (weight-GETWT(tail,head))*m->workspace[i];
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
