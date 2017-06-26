#include "ergm_changestat_multinet.h"
#include "ergm_changestat.h"
#include "ergm_changestat_operator.h"


I_CHANGESTAT_FN(i__subnets){
  double *inputs = INPUT_PARAM;
  ALLOC_AUX_STORAGE(1, StoreSubnets, sn); inputs++;
  sn->ns = *(inputs++);
  sn->inwp = nwp;
  sn->onwp = Calloc(sn->ns, Network);

  /* Set up the layer information. */
  sn->sid = inputs - 1; // The -1 is because Vertex IDs count from 1.
  inputs += N_NODES;
  sn->smap = inputs - 1;
  inputs += N_NODES;

  for(unsigned int i=0; i<sn->ns; i++){
    Vertex lnnodes, lbip;
    if(BIPARTITE){
      lbip = lnnodes = *(inputs++);
      lnnodes += *(inputs++);
    }else{
      lbip = 0;
      lnnodes = *(inputs++);
    }

    sn->onwp[i] = NetworkInitialize(NULL, NULL, 0, lnnodes, DIRECTED, lbip, 0, 0, NULL);
  }
  
  EXEC_THROUGH_NET_EDGES(t, h, e, {
      ToggleEdge(MN_IO_TAIL(sn, t), MN_IO_HEAD(sn, h), sn->onwp + MN_SID_TAIL(sn, t));
    });
}

U_CHANGESTAT_FN(u__subnets){ 
  GET_AUX_STORAGE(StoreSubnets, sn);
  ToggleEdge(MN_IO_TAIL(sn, tail), MN_IO_HEAD(sn, head),sn->onwp + MN_SID_TAIL(sn, tail));
}

F_CHANGESTAT_FN(f__subnets){
  GET_AUX_STORAGE(StoreSubnets, sn);
  for(unsigned int i=0; i<sn->ns; i++)
    NetworkDestroy(sn->onwp + i);
  Free(sn->onwp);
}

I_CHANGESTAT_FN(i_MultiNet){
  double *inputs = INPUT_PARAM; 
  GET_AUX_STORAGE(StoreSubnets, sn); inputs++;
  unsigned int ns = sn->ns;
  double *w = inputs;
  inputs+=ns;
  ALLOC_STORAGE(ns, Model*, ms);

  for(unsigned int i=0; i<sn->ns; i++){
    ms[i] = unpack_Model_as_double(&inputs);
    InitStats(sn->onwp + i, ms[i]);
  }
}

C_CHANGESTAT_FN(c_MultiNet){
  GET_AUX_STORAGE(StoreSubnets, sn);
  GET_STORAGE(Model*, ms);
  double *w = INPUT_PARAM+1;

  for(unsigned int i=0; i<sn->ns; i++){
    if(w[i]){
      Vertex st = MN_IO_TAIL(sn, tail), sh = MN_IO_HEAD(sn, head);
      ChangeStats(1, &st, &sh, sn->onwp + i, ms[i]);
      for(unsigned int j=0; j<N_CHANGE_STATS; j++)
	CHANGE_STAT[j] += ms[i]->workspace[j]*w[i];
    }
  }
}

U_CHANGESTAT_FN(u_MultiNet){
  GET_AUX_STORAGE(StoreSubnets, sn);
  GET_STORAGE(Model*, ms);
  double *w = INPUT_PARAM+1;

  for(unsigned int i=0; i<sn->ns; i++){
    Vertex st = MN_IO_TAIL(sn, tail), sh = MN_IO_HEAD(sn, head);
    Model *m = ms[i];
    UPDATE_STORAGE(st, sh, sn->onwp + i, m, NULL);
  }
}

F_CHANGESTAT_FN(f_MultiNet){
  GET_AUX_STORAGE(StoreSubnets, sn);
  GET_STORAGE(Model*, ms);
  double *w = INPUT_PARAM+1;

  for(unsigned int i=0; i<sn->ns; i++){
    ModelDestroy(sn->onwp + i, ms[i]);
  }
}

