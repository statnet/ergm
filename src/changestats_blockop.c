#include "ergm_changestats_blockop.h"
#include "ergm_changestats_auxnet.h"

/* NodematchFilter(formula) */

I_CHANGESTAT_FN(i_NodematchFilter){
  GET_AUX_STORAGE(StoreNetAndBID, blkinfo);
  double *inputs = INPUT_PARAM+1;
  Network *bnwp=blkinfo->nwp;
  // No need to allocate it: we are only storing a pointer to a model.
  STORAGE = unpack_Model_as_double(&inputs);
 
  InitStats(bnwp, STORAGE);
}

C_CHANGESTAT_FN(c_NodematchFilter){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(StoreNetAndBID, blkinfo);
  Network *bnwp=blkinfo->nwp;
  double *b = blkinfo->b;

  if(b[tail]==b[head]){
    ChangeStats(1, &tail, &head, bnwp, m);
    memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double));
  }
}

U_CHANGESTAT_FN(u_NodematchFilter){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(StoreNetAndBID, blkinfo);
  Network *bnwp=blkinfo->nwp;
  double *b = blkinfo->b;

  if(b[tail]==b[head])
    UPDATE_STORAGE(tail, head, bnwp, m, NULL);
}

F_CHANGESTAT_FN(f_NodematchFilter){
  GET_STORAGE(Model, m);

  ModelDestroy(nwp, m);

  STORAGE = NULL;
}

