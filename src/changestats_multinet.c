#include "ergm_changestat_multinet.h"
#include "ergm_changestat.h"
#include "ergm_changestat_operator.h"


I_CHANGESTAT_FN(i__subnets){
  double *inputs = INPUT_PARAM;
  ALLOC_AUX_STORAGE(1, StoreSubnets, sn); inputs++;
  sn->ns = *(inputs++);
  sn->inwp = nwp;
  sn->onwp = Calloc(sn->ns, Network);
  sn->onwp--; // The -- is because Network IDs count from 1.

  /* Set up the layer information. */
  sn->sid = inputs - 1; // The -1 is because Vertex IDs count from 1.
  inputs += N_NODES;
  sn->smap = inputs - 1;
  inputs += N_NODES;

  for(unsigned int i=1; i<=sn->ns; i++){
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
  
  EXEC_THROUGH_NET_EDGES_PRE(t, h, e, {
      ToggleEdge(MN_IO_TAIL(sn, t), MN_IO_HEAD(sn, h), sn->onwp + MN_SID_TAIL(sn, t));
    });
}

U_CHANGESTAT_FN(u__subnets){ 
  GET_AUX_STORAGE(StoreSubnets, sn);
  ToggleEdge(MN_IO_TAIL(sn, tail), MN_IO_HEAD(sn, head),sn->onwp + MN_SID_TAIL(sn, tail));
}

F_CHANGESTAT_FN(f__subnets){
  GET_AUX_STORAGE(StoreSubnets, sn);
  for(unsigned int i=1; i<=sn->ns; i++)
    NetworkDestroy(sn->onwp + i);
  sn->onwp++;
  Free(sn->onwp);
}

// MultiNet: Take a weighted networkwise sum of the networks' statistics.

I_CHANGESTAT_FN(i_MultiNet){
  double *inputs = INPUT_PARAM; 
  GET_AUX_STORAGE(StoreSubnets, sn); inputs++;
  unsigned int ns = sn->ns;
  inputs+=ns;
  ALLOC_STORAGE(ns, Model*, ms);

  for(unsigned int i=1; i<=sn->ns; i++){
    ms[i-1] = unpack_Model_as_double(&inputs);
    InitStats(sn->onwp + i, ms[i-1]);
  }
}

C_CHANGESTAT_FN(c_MultiNet){
  GET_AUX_STORAGE(StoreSubnets, sn);
  GET_STORAGE(Model*, ms);
  double *w = INPUT_PARAM; // Subnetworks are coded from 1.

  unsigned int i = MN_SID_TAIL(sn, tail);
  if(w[i]){
    Vertex st = MN_IO_TAIL(sn, tail), sh = MN_IO_HEAD(sn, head);
    Model *m = ms[i-1];
      ChangeStats(1, &st, &sh, sn->onwp + i, m);
      for(unsigned int j=0; j<N_CHANGE_STATS; j++)
	CHANGE_STAT[j] += m->workspace[j]*w[i];
  }
}

U_CHANGESTAT_FN(u_MultiNet){
  GET_AUX_STORAGE(StoreSubnets, sn);
  GET_STORAGE(Model*, ms);

  unsigned int i = MN_SID_TAIL(sn, tail);
  Vertex st = MN_IO_TAIL(sn, tail), sh = MN_IO_HEAD(sn, head);
  UPDATE_STORAGE(st, sh, sn->onwp + i, ms[i-1], NULL);
}

F_CHANGESTAT_FN(f_MultiNet){
  GET_AUX_STORAGE(StoreSubnets, sn);
  GET_STORAGE(Model*, ms);

  for(unsigned int i=1; i<=sn->ns; i++){
    ModelDestroy(sn->onwp + i, ms[i-1]);
  }
}

// MultiNets: Concatenate the networks' statistics; network statistic counts may be heterogeneous.

I_CHANGESTAT_FN(i_MultiNets){
  double *inputs = INPUT_PARAM; 
  GET_AUX_STORAGE(StoreSubnets, sn); inputs++;
  unsigned int ns = sn->ns;
  inputs+=ns;
  ALLOC_STORAGE(ns, Model*, ms);

  for(unsigned int i=1; i<=sn->ns; i++){
    ms[i-1] = unpack_Model_as_double(&inputs);
    InitStats(sn->onwp + i, ms[i-1]);
  }
}

C_CHANGESTAT_FN(c_MultiNets){
  double *pos = INPUT_PARAM; // Starting positions of subnetworks' statistics.
  GET_AUX_STORAGE(StoreSubnets, sn); pos++;
  GET_STORAGE(Model*, ms);

  unsigned int i = MN_SID_TAIL(sn, tail);
  Vertex st = MN_IO_TAIL(sn, tail), sh = MN_IO_HEAD(sn, head);
  Model *m = ms[i-1];
  ChangeStats(1, &st, &sh, sn->onwp + i, m);
  memcpy(CHANGE_STAT + (unsigned int)(pos[i-1]), m->workspace, m->n_stats*sizeof(double));
}

U_CHANGESTAT_FN(u_MultiNets){
  GET_AUX_STORAGE(StoreSubnets, sn);
  GET_STORAGE(Model*, ms);

  unsigned int i = MN_SID_TAIL(sn, tail);
  Vertex st = MN_IO_TAIL(sn, tail), sh = MN_IO_HEAD(sn, head);
  UPDATE_STORAGE(st, sh, sn->onwp + i, ms[i-1], NULL);
}

F_CHANGESTAT_FN(f_MultiNets){
  GET_AUX_STORAGE(StoreSubnets, sn);
  GET_STORAGE(Model*, ms);

  for(unsigned int i=1; i<=sn->ns; i++){
    ModelDestroy(sn->onwp + i, ms[i-1]);
  }
}

