#include "changestat_operator.h"

typedef struct{StoreAuxAndModel am; unsigned int n_stats_1, n_stats_2;} StoreAuxAndModelAnd2Stats;

I_CHANGESTAT_FN(i_passthrough_term){
  double *inputs = INPUT_PARAM;
  ALLOC_STORAGE(1, StoreAuxAndModel, store);

  Model *m = store->m = unpack_Modelasdouble(&inputs);

  // Before we pass nwp down, we swap the aux_storage in for our own version.
  void **nwp_aux_storage = nwp->aux_storage;
  if(m->n_aux){
    store->aux_storage = nwp->aux_storage = (void **)malloc(sizeof(void *)*m->n_aux);
    for(unsigned int i = 0; i<m->n_aux; i++) nwp->aux_storage[i] = NULL;
  }else nwp->aux_storage = NULL;

  InitStats(nwp, m);

  // Now put nwp's aux_storage back:
  nwp->aux_storage = nwp_aux_storage;
}

D_CHANGESTAT_FN(d_passthrough_term){
  GET_STORAGE(StoreAuxAndModel, store);
  Model *m = store->m;
    
  void **nwp_aux_storage = nwp->aux_storage;
  nwp->aux_storage = store->aux_storage;
  ChangeStats(ntoggles, tails, heads, nwp, m);
  nwp->aux_storage = nwp_aux_storage;
  
  memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double));
}

U_CHANGESTAT_FN(u_passthrough_term){
  GET_STORAGE(StoreAuxAndModel, store);
  Model *m = store->m;

  void **nwp_aux_storage = nwp->aux_storage;
  nwp->aux_storage = store->aux_storage;
  UPDATE_STORAGE(tail, head, m, nwp);
  nwp->aux_storage = nwp_aux_storage;
}

F_CHANGESTAT_FN(f_passthrough_term){
  GET_STORAGE(StoreAuxAndModel, store);
  Model *m = store->m;

  void **nwp_aux_storage = nwp->aux_storage;
  nwp->aux_storage = store->aux_storage;
  ModelDestroy(m, nwp);
  nwp->aux_storage = nwp_aux_storage;
}


I_CHANGESTAT_FN(i_interact){
  double *inputs = INPUT_PARAM;
  ALLOC_STORAGE(1, StoreAuxAndModelAnd2Stats, store);

  store->n_stats_1 = *(inputs++);
  store->n_stats_2 = *(inputs++);
  Model *m = store->am.m = unpack_Modelasdouble(&inputs);

  // Before we pass nwp down, we swap the aux_storage in for our own version.
  void **nwp_aux_storage = nwp->aux_storage;
  if(m->n_aux){
    store->am.aux_storage = nwp->aux_storage = (void **)malloc(sizeof(void *)*m->n_aux);
    for(unsigned int i = 0; i<m->n_aux; i++) nwp->aux_storage[i] = NULL;
  }else nwp->aux_storage = NULL;

  InitStats(nwp, m);

  // Now put nwp's aux_storage back:
  nwp->aux_storage = nwp_aux_storage;
}

C_CHANGESTAT_FN(c_interact){
  
  GET_STORAGE(StoreAuxAndModelAnd2Stats, store);
  Model *m = store->am.m;
  
  void **nwp_aux_storage = nwp->aux_storage;
  nwp->aux_storage = store->am.aux_storage;
  ChangeStats(1, &tail, &head, nwp, m);
  nwp->aux_storage = nwp_aux_storage;

  double *w2 = m->workspace + store->n_stats_1;
  unsigned int pos = 0;
  for(unsigned int j=0; j<store->n_stats_2; j++){
    for(unsigned int i=0; i<store->n_stats_1; i++){
      CHANGE_STAT[pos++] = m->workspace[i]*w2[j];
    }
  }
}

U_CHANGESTAT_FN(u_interact){
  GET_STORAGE(StoreAuxAndModelAnd2Stats, store);
  Model *m = store->am.m;

  void **nwp_aux_storage = nwp->aux_storage;
  nwp->aux_storage = store->am.aux_storage;
  UPDATE_STORAGE(tail, head, m, nwp);
  nwp->aux_storage = nwp_aux_storage;
}

F_CHANGESTAT_FN(f_interact){
  GET_STORAGE(StoreAuxAndModelAnd2Stats, store);
  Model *m = store->am.m;

  void **nwp_aux_storage = nwp->aux_storage;
  nwp->aux_storage = store->am.aux_storage;
  ModelDestroy(m, nwp);
  nwp->aux_storage = nwp_aux_storage;
}


I_CHANGESTAT_FN(i_main_interact){
  double *inputs = INPUT_PARAM;
  ALLOC_STORAGE(1, StoreAuxAndModelAnd2Stats, store);

  store->n_stats_1 = *(inputs++);
  store->n_stats_2 = *(inputs++);
  Model *m = store->am.m = unpack_Modelasdouble(&inputs);

  // Before we pass nwp down, we swap the aux_storage in for our own version.
  void **nwp_aux_storage = nwp->aux_storage;
  if(m->n_aux){
    store->am.aux_storage = nwp->aux_storage = (void **)malloc(sizeof(void *)*m->n_aux);
    for(unsigned int i = 0; i<m->n_aux; i++) nwp->aux_storage[i] = NULL;
  }else nwp->aux_storage = NULL;

  InitStats(nwp, m);

  // Now put nwp's aux_storage back:
  nwp->aux_storage = nwp_aux_storage;
}

C_CHANGESTAT_FN(c_main_interact){
  
  GET_STORAGE(StoreAuxAndModelAnd2Stats, store);
  Model *m = store->am.m;
  
  void **nwp_aux_storage = nwp->aux_storage;
  nwp->aux_storage = store->am.aux_storage;
  ChangeStats(1, &tail, &head, nwp, m);
  nwp->aux_storage = nwp_aux_storage;

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
  GET_STORAGE(StoreAuxAndModelAnd2Stats, store);
  Model *m = store->am.m;

  void **nwp_aux_storage = nwp->aux_storage;
  nwp->aux_storage = store->am.aux_storage;
  UPDATE_STORAGE(tail, head, m, nwp);
  nwp->aux_storage = nwp_aux_storage;
}

F_CHANGESTAT_FN(f_main_interact){
  GET_STORAGE(StoreAuxAndModelAnd2Stats, store);
  Model *m = store->am.m;

  void **nwp_aux_storage = nwp->aux_storage;
  nwp->aux_storage = store->am.aux_storage;
  ModelDestroy(m, nwp);
  nwp->aux_storage = nwp_aux_storage;
}
