#include "ergm_model.h"
#include "ergm_storage.h"
#include "ergm_changestat_operator.h"

typedef struct{Model *m; unsigned int n_stats_1, n_stats_2;} StoreModelAnd2Stats;

I_CHANGESTAT_FN(i_interact){
  double *inputs = INPUT_PARAM;
  ALLOC_STORAGE(1, StoreModelAnd2Stats, store);

  store->n_stats_1 = *(inputs++);
  store->n_stats_2 = *(inputs++);
  store->m = ModelInitialize(getListElement(mtp->R, "submodel"), NULL,  nwp, FALSE);
  DELETE_IF_UNUSED_IN_SUBMODEL(u_func, store->m);
  DELETE_IF_UNUSED_IN_SUBMODEL(z_func, store->m);
}

C_CHANGESTAT_FN(c_interact){
  GET_STORAGE(StoreModelAnd2Stats, store);
  Model *m = store->m;

  int change = edgeflag ? -1 : +1;
  ChangeStats1(tail, head, nwp, m, edgeflag);

  double *w2 = m->workspace + store->n_stats_1;
  unsigned int pos = 0;
  for(unsigned int j=0; j<store->n_stats_2; j++){
    for(unsigned int i=0; i<store->n_stats_1; i++){
      CHANGE_STAT[pos++] = m->workspace[i]*w2[j]*change;
    }
  }
}

U_CHANGESTAT_FN(u_interact){
  GET_STORAGE(StoreModelAnd2Stats, store);
  Model *m = store->m;

  UPDATE_STORAGE(tail, head, nwp, m, NULL, edgeflag);
}

Z_CHANGESTAT_FN(z_interact){
  GET_STORAGE(StoreModelAnd2Stats, store);
  Model *m = store->m;

  ZStats(nwp, m);

  double *w2 = m->workspace + store->n_stats_1;
  unsigned int pos = 0;
  for(unsigned int j=0; j<store->n_stats_2; j++){
    for(unsigned int i=0; i<store->n_stats_1; i++){
      CHANGE_STAT[pos++] = m->workspace[i]*w2[j];
    }
  }
}

F_CHANGESTAT_FN(f_interact){
  GET_STORAGE(StoreModelAnd2Stats, store);
  Model *m = store->m;

  ModelDestroy(nwp, m);
}


I_CHANGESTAT_FN(i_main_interact){
  double *inputs = INPUT_PARAM;
  ALLOC_STORAGE(1, StoreModelAnd2Stats, store);

  store->n_stats_1 = *(inputs++);
  store->n_stats_2 = *(inputs++);
  store->m = ModelInitialize(getListElement(mtp->R, "submodel"), NULL,  nwp, FALSE);
  DELETE_IF_UNUSED_IN_SUBMODEL(u_func, store->m);
  DELETE_IF_UNUSED_IN_SUBMODEL(z_func, store->m);
}

C_CHANGESTAT_FN(c_main_interact){
  GET_STORAGE(StoreModelAnd2Stats, store);
  Model *m = store->m;
  
  int change = edgeflag ? -1 : +1;
  ChangeStats1(tail, head, nwp, m, edgeflag);

  double *w2 = m->workspace + store->n_stats_1;
  unsigned int pos = 0;

  for(unsigned int i=0; i<store->n_stats_1; i++){
    CHANGE_STAT[pos++] = m->workspace[i];
  }
  for(unsigned int j=0; j<store->n_stats_2; j++){
    CHANGE_STAT[pos++] = w2[j];
  }
  for(unsigned int j=0; j<store->n_stats_2; j++){
    for(unsigned int i=0; i<store->n_stats_1; i++){
      CHANGE_STAT[pos++] = m->workspace[i]*w2[j] * change;
    }
  }
}

U_CHANGESTAT_FN(u_main_interact){
  GET_STORAGE(StoreModelAnd2Stats, store);
  Model *m = store->m;

  UPDATE_STORAGE(tail, head, nwp, m, NULL, edgeflag);
}

Z_CHANGESTAT_FN(z_main_interact){
  GET_STORAGE(StoreModelAnd2Stats, store);
  Model *m = store->m;

  ZStats(nwp, m);

  double *w2 = m->workspace + store->n_stats_1;
  unsigned int pos = 0;

  for(unsigned int i=0; i<store->n_stats_1; i++){
    CHANGE_STAT[pos++] = m->workspace[i];
  }
  for(unsigned int j=0; j<store->n_stats_2; j++){
    CHANGE_STAT[pos++] = w2[j];
  }
  for(unsigned int j=0; j<store->n_stats_2; j++){
    for(unsigned int i=0; i<store->n_stats_1; i++){
      CHANGE_STAT[pos++] = m->workspace[i]*w2[j];
    }
  }
}

F_CHANGESTAT_FN(f_main_interact){
  GET_STORAGE(StoreModelAnd2Stats, store);
  Model *m = store->m;

  ModelDestroy(nwp, m);

  STORAGE=NULL;
}
