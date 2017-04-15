#include "changestats_blockop.h"

/* within_block(formula) */

I_CHANGESTAT_FN(i_within_block){
  GET_AUX_STORAGE(StoreNetAndBID, blkinfo);
  double *inputs = INPUT_PARAM+1;
  Network *bnwp=&(blkinfo->nw);
  // No need to allocate it: we are only storing a pointer to a model.
  STORAGE = unpack_Model_as_double(&inputs);
 
  InitStats(bnwp, STORAGE);
}

C_CHANGESTAT_FN(c_within_block){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(StoreNetAndBID, blkinfo);
  Network *bnwp=&(blkinfo->nw);
  double *b = blkinfo->b;

  if(b[tail]==b[head]){
    ChangeStats(1, &tail, &head, bnwp, m);
    memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double));
  }
}

U_CHANGESTAT_FN(u_within_block){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(StoreNetAndBID, blkinfo);
  Network *bnwp=&(blkinfo->nw);
  double *b = blkinfo->b;

  if(b[tail]==b[head])
    UPDATE_STORAGE(tail, head, bnwp, m, NULL);
}

F_CHANGESTAT_FN(f_within_block){
  GET_STORAGE(Model, m);

  ModelDestroy(nwp, m);

  STORAGE = NULL;
}

/* within_block:

   maintains a network that mirrors the main network but excludes ties
   not within specified blocks; also exports a numeric vector of block
   memberships
 */

I_CHANGESTAT_FN(i__within_block){
  ALLOC_AUX_STORAGE(1, StoreNetAndBID, blkinfo);
  blkinfo->nw = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, 0, 0, NULL);
  Network *bnwp = &(blkinfo->nw);
  double *b = blkinfo->b = INPUT_PARAM;

  for(Vertex tail=1; tail <= N_TAILS; tail++){
    Vertex head;
    Edge e;
    STEP_THROUGH_OUTEDGES(tail, e, head) {
      if(b[tail]==b[head])
	ToggleEdge(tail, head, bnwp);
    }
  }
}

U_CHANGESTAT_FN(u__within_block){
  GET_AUX_STORAGE(StoreNetAndBID, blkinfo);
  Network *bnwp = &(blkinfo->nw);
  double *b = blkinfo->b;

  if(b[tail]==b[head])
    ToggleEdge(tail, head, bnwp);
}

F_CHANGESTAT_FN(f__within_block){
  GET_AUX_STORAGE(StoreNetAndBID, blkinfo);
  Network *bnwp = &(blkinfo->nw);
  NetworkDestroy(bnwp);
}
