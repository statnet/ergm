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

