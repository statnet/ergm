#include "changestat.h"
#include "model.h"

typedef struct{void **aux_storage; Model *m;} StoreAuxAndModel;

static inline unsigned char *dec_str(double **x){
  unsigned int l = (*x)[0];
  unsigned char *s = (unsigned char *) malloc((l+1)*sizeof(char));
  for(unsigned int i=0; i<l; i++){
    s[i] = (unsigned char) (unsigned int) (*x)[i+1];
  }
  s[l] = (unsigned char) 0;
  (*x)+=l+1;
  return s;
}

I_CHANGESTAT_FN(i_meta_term){
  double *inputs = INPUT_PARAM;
  ALLOC_STORAGE(1, StoreAuxAndModel, store);

  int n_terms = *(inputs++);
  char *fnames = (char *) dec_str(&inputs);
  char *snames = (char *) dec_str(&inputs);
  
  Model *m = store->m = ModelInitialize(fnames, snames, &inputs, n_terms);
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

D_CHANGESTAT_FN(d_meta_term){
  GET_STORAGE(StoreAuxAndModel, store);
  Model *m = store->m;
    
  void **nwp_aux_storage = nwp->aux_storage;
  nwp->aux_storage = store->aux_storage;
  ChangeStats(ntoggles, tails, heads, nwp, m);
  nwp->aux_storage = nwp_aux_storage;
  
  memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double));
}

U_CHANGESTAT_FN(u_meta_term){
  GET_STORAGE(StoreAuxAndModel, store);
  Model *m = store->m;

  void **nwp_aux_storage = nwp->aux_storage;
  nwp->aux_storage = store->aux_storage;
  UPDATE_STORAGE(tail, head, m, nwp);
  nwp->aux_storage = nwp_aux_storage;
}

F_CHANGESTAT_FN(f_meta_term){
  GET_STORAGE(StoreAuxAndModel, store);
  Model *m = store->m;

  void **nwp_aux_storage = nwp->aux_storage;
  nwp->aux_storage = store->aux_storage;
  ModelDestroy(m, nwp);
  nwp->aux_storage = nwp_aux_storage;
}


I_CHANGESTAT_FN(i_interact){
  double *inputs = INPUT_PARAM;
  ALLOC_STORAGE(1, StoreAuxAndModel, store);

  int n_terms = 2;
  char *fnames = (char *) dec_str(&inputs);
  char *snames = (char *) dec_str(&inputs);
  
  Model *m = store->m = ModelInitialize(fnames, snames, &inputs, n_terms);
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

C_CHANGESTAT_FN(c_interact){
  GET_STORAGE(StoreAuxAndModel, store);
  Model *m = store->m;
    
  void **nwp_aux_storage = nwp->aux_storage;
  nwp->aux_storage = store->aux_storage;
  ChangeStats(1, &tail, &head, nwp, m);
  nwp->aux_storage = nwp_aux_storage;

  double *w2 = m->workspace + m->termarray[0].nstats;
  unsigned int pos = 0;
  for(unsigned int j=0; j<m->termarray[1].nstats; j++){
    for(unsigned int i=0; i<m->termarray[0].nstats; i++){
      CHANGE_STAT[pos++] = m->workspace[i]*w2[j];
    }
  }
}

U_CHANGESTAT_FN(u_interact){
  GET_STORAGE(StoreAuxAndModel, store);
  Model *m = store->m;

  void **nwp_aux_storage = nwp->aux_storage;
  nwp->aux_storage = store->aux_storage;
  UPDATE_STORAGE(tail, head, m, nwp);
  nwp->aux_storage = nwp_aux_storage;
}

F_CHANGESTAT_FN(f_interact){
  GET_STORAGE(StoreAuxAndModel, store);
  Model *m = store->m;

  void **nwp_aux_storage = nwp->aux_storage;
  nwp->aux_storage = store->aux_storage;
  ModelDestroy(m, nwp);
  nwp->aux_storage = nwp_aux_storage;
}


I_CHANGESTAT_FN(i_main_interact){
  double *inputs = INPUT_PARAM;
  ALLOC_STORAGE(1, StoreAuxAndModel, store);

  int n_terms = 2;
  char *fnames = (char *) dec_str(&inputs);
  char *snames = (char *) dec_str(&inputs);
  
  Model *m = store->m = ModelInitialize(fnames, snames, &inputs, n_terms);
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

C_CHANGESTAT_FN(c_main_interact){
  GET_STORAGE(StoreAuxAndModel, store);
  Model *m = store->m;
    
  void **nwp_aux_storage = nwp->aux_storage;
  nwp->aux_storage = store->aux_storage;
  ChangeStats(1, &tail, &head, nwp, m);
  nwp->aux_storage = nwp_aux_storage;

  double *w2 = m->workspace + m->termarray[0].nstats;
  unsigned int pos = 0;

  for(unsigned int i=0; i<m->termarray[0].nstats; i++){
    CHANGE_STAT[pos++] = m->workspace[i];
  }
  for(unsigned int j=0; j<m->termarray[1].nstats; j++){
    CHANGE_STAT[pos++] = w2[j];
  }
  for(unsigned int j=0; j<m->termarray[1].nstats; j++){
    for(unsigned int i=0; i<m->termarray[0].nstats; i++){
      CHANGE_STAT[pos++] = m->workspace[i]*w2[j];
    }
  }
}

U_CHANGESTAT_FN(u_main_interact){
  GET_STORAGE(StoreAuxAndModel, store);
  Model *m = store->m;

  void **nwp_aux_storage = nwp->aux_storage;
  nwp->aux_storage = store->aux_storage;
  UPDATE_STORAGE(tail, head, m, nwp);
  nwp->aux_storage = nwp_aux_storage;
}

F_CHANGESTAT_FN(f_main_interact){
  GET_STORAGE(StoreAuxAndModel, store);
  Model *m = store->m;

  void **nwp_aux_storage = nwp->aux_storage;
  nwp->aux_storage = store->aux_storage;
  ModelDestroy(m, nwp);
  nwp->aux_storage = nwp_aux_storage;
}
