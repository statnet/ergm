#include "ergm_changestat_operator.h"
#include "ergm_changestat_auxnet.h"
#include "ergm_util.h"

Model *unpack_Model_as_double(double **x){
  int n_terms = *((*x)++);
  char *fnames = (char *) unpack_str_as_double(x);
  char *snames = (char *) unpack_str_as_double(x);
  Model *m = ModelInitialize(fnames, snames, x, n_terms);
  Free(fnames);
  Free(snames);
  return m;  
}

/* OnAuxnet: A generic implementation that evaluates a model on an
   auxiliary network. Covers many special cases. */
I_CHANGESTAT_FN(i_OnAuxnet){
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  double *inputs = INPUT_PARAM + 1;
  Model *m = STORAGE = unpack_Model_as_double(&inputs);
  InitStats(auxnet->onwp, m);
}

C_CHANGESTAT_FN(c_OnAuxnet){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(StoreAuxnet, auxnet);

  MAP_TOGGLE_1(tail, head, auxnet, ntoggles, tails, heads);

  if(ntoggles){
    ChangeStats(1, tails, heads, auxnet->onwp, m);
    memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double));
  }
}

U_CHANGESTAT_FN(u_OnAuxnet){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(StoreAuxnet, auxnet);

  MAP_TOGGLE_1(tail, head, auxnet, ntoggles, tails, heads);
  if(ntoggles) UPDATE_STORAGE(*tails, *heads, auxnet->onwp, m, NULL, edgeflag);
}

F_CHANGESTAT_FN(f_OnAuxnet){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(StoreAuxnet, auxnet);
  ModelDestroy(auxnet->onwp, m);
  STORAGE = NULL;
}
