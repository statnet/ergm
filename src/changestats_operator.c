#include "changestat_operator.h"

typedef struct{Model *m; unsigned int n_stats_1, n_stats_2;} StoreModelAnd2Stats;

I_CHANGESTAT_FN(i_passthrough_term){
  double *inputs = INPUT_PARAM;
  GET_STORAGE(Model, m); // No need to allocate it: we are only storing a pointer to a model.

  m = unpack_Modelasdouble(&inputs);

  InitStats(nwp, m);
}

D_CHANGESTAT_FN(d_passthrough_term){
  GET_STORAGE(Model, m);

  ChangeStats(ntoggles, tails, heads, nwp, m);
  
  memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double));
}

U_CHANGESTAT_FN(u_passthrough_term){
  GET_STORAGE(Model, m);

  UPDATE_STORAGE(tail, head, m, nwp);
}

F_CHANGESTAT_FN(f_passthrough_term){
  GET_STORAGE(Model, m);

  ModelDestroy(m, nwp);

  mtp->storage=NULL;
}


I_CHANGESTAT_FN(i_interact){
  double *inputs = INPUT_PARAM;
  ALLOC_STORAGE(1, StoreModelAnd2Stats, store);

  store->n_stats_1 = *(inputs++);
  store->n_stats_2 = *(inputs++);
  Model *m = store->m = unpack_Modelasdouble(&inputs);

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

  UPDATE_STORAGE(tail, head, m, nwp);
}

F_CHANGESTAT_FN(f_interact){
  GET_STORAGE(StoreModelAnd2Stats, store);
  Model *m = store->m;

  ModelDestroy(m, nwp);
}


I_CHANGESTAT_FN(i_main_interact){
  double *inputs = INPUT_PARAM;
  ALLOC_STORAGE(1, StoreModelAnd2Stats, store);

  store->n_stats_1 = *(inputs++);
  store->n_stats_2 = *(inputs++);
  Model *m = store->m = unpack_Modelasdouble(&inputs);

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

  UPDATE_STORAGE(tail, head, m, nwp);
}

F_CHANGESTAT_FN(f_main_interact){
  GET_STORAGE(StoreModelAnd2Stats, store);
  Model *m = store->m;

  ModelDestroy(m, nwp);

  mtp->storage=NULL;
}
