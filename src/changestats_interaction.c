#include "ergm_changestat_operator.h"

typedef struct{Model *m; unsigned int n_stats_1, n_stats_2;} StoreModelAnd2Stats;

I_CHANGESTAT_FN(i_interact){
  double *inputs = INPUT_PARAM;
  ALLOC_STORAGE(1, StoreModelAnd2Stats, store);

  store->n_stats_1 = *(inputs++);
  store->n_stats_2 = *(inputs++);
  Model *m = store->m = unpack_Model_as_double(&inputs);

  InitStats(nwp, m);
}

C_CHANGESTAT_FN(c_interact){
  GET_STORAGE(StoreModelAnd2Stats, store);
  Model *m = store->m;
  
  ChangeStats(1, &tail, &head, nwp, m);

  double *w2 = m->workspace + store->n_stats_1;
  unsigned int pos = 0;
  for(unsigned int j=0; j<store->n_stats_2; j++){
    for(unsigned int i=0; i<store->n_stats_1; i++){
      CHANGE_STAT[pos++] = m->workspace[i]*w2[j];
    }
  }
}

U_CHANGESTAT_FN(u_interact){
  GET_STORAGE(StoreModelAnd2Stats, store);
  Model *m = store->m;

  UPDATE_STORAGE(tail, head, nwp, m, NULL);
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
  Model *m = store->m = unpack_Model_as_double(&inputs);

  InitStats(nwp, m);
}

C_CHANGESTAT_FN(c_main_interact){
  
  GET_STORAGE(StoreModelAnd2Stats, store);
  Model *m = store->m;
  
  ChangeStats(1, &tail, &head, nwp, m);

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

U_CHANGESTAT_FN(u_main_interact){
  GET_STORAGE(StoreModelAnd2Stats, store);
  Model *m = store->m;

  UPDATE_STORAGE(tail, head, nwp, m, NULL);
}

F_CHANGESTAT_FN(f_main_interact){
  GET_STORAGE(StoreModelAnd2Stats, store);
  Model *m = store->m;

  ModelDestroy(nwp, m);

  STORAGE=NULL;
}
