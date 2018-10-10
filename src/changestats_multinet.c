#include "ergm_changestat_multinet.h"
#include "ergm_changestat.h"
#include "ergm_changestat_operator.h"


I_CHANGESTAT_FN(i__subnets){
  double *inputs = INPUT_PARAM;
  ALLOC_AUX_STORAGE(1, StoreSubnets, sn); inputs++;
  sn->ns = *(inputs++);
  sn->inwp = nwp;
  sn->onwp = Calloc(sn->ns, Network *);
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
      ToggleEdge(MN_IO_TAIL(sn, t), MN_IO_HEAD(sn, h), sn->onwp[MN_SID_TAIL(sn, t)]);
    });
}

U_CHANGESTAT_FN(u__subnets){ 
  GET_AUX_STORAGE(StoreSubnets, sn);
  ToggleEdge(MN_IO_TAIL(sn, tail), MN_IO_HEAD(sn, head),sn->onwp[MN_SID_TAIL(sn, tail)]);
}

F_CHANGESTAT_FN(f__subnets){
  GET_AUX_STORAGE(StoreSubnets, sn);
  for(unsigned int i=1; i<=sn->ns; i++)
    NetworkDestroy(sn->onwp[i]);
  sn->onwp++;
  Free(sn->onwp);
}

// MultiNet: Take a weighted networkwise sum of the networks' statistics.

I_CHANGESTAT_FN(i_MultiNet){
  /*
    inputs expects:
    1: position of the subnets auxiliary
    1: number of weights per network (nwts)
    nwts*ns: matrix of weights, in network-major order
    ?*ns: submodel specifications for nm submodels
  */
  
  double *inputs = INPUT_PARAM; 
  GET_AUX_STORAGE(StoreSubnets, sn); inputs++;
  unsigned int ns = sn->ns;
  unsigned int nwts = *(inputs++);
  double *wts = inputs; inputs+=ns*nwts;
  
  ALLOC_STORAGE(ns, Model*, ms);

  for(unsigned int i=1; i<=sn->ns; i++){
    unsigned int used=FALSE;
    for(unsigned int j=0; j<nwts; j++){
      if(wts[j]!=0){
	used=TRUE;
	break;
      }
    }
    wts += nwts; // OK to clobber it here.
    if(used){
      ms[i-1] = unpack_Model_as_double(&inputs);
      InitStats(sn->onwp[i], ms[i-1]);
    }else ms[i-1] = NULL;
  }
}

C_CHANGESTAT_FN(c_MultiNet){
  double *inputs = INPUT_PARAM; 
  GET_AUX_STORAGE(StoreSubnets, sn); inputs++;
  GET_STORAGE(Model*, ms);
  unsigned int nwts = *(inputs++);
  double *wts = inputs;

  unsigned int i = MN_SID_TAIL(sn, tail);
  Model *m = ms[i-1];
  if(m){ // NULL if network has weights 0.
    Vertex st = MN_IO_TAIL(sn, tail), sh = MN_IO_HEAD(sn, head);
    ChangeStats(1, &st, &sh, sn->onwp[i], m);

    wts += (i-1)*nwts; // Position of that network's weight vector.
    for(unsigned int j=0; j<m->n_stats; j++)
      for(unsigned int k=0; k<nwts; k++)
	CHANGE_STAT[j*nwts+k] += m->workspace[j]*wts[k];
  }
}

U_CHANGESTAT_FN(u_MultiNet){
  GET_AUX_STORAGE(StoreSubnets, sn);
  GET_STORAGE(Model*, ms);

  unsigned int i = MN_SID_TAIL(sn, tail);
  Model *m = ms[i-1];
  if(m){ // NULL if network has weights 0.
    Vertex st = MN_IO_TAIL(sn, tail), sh = MN_IO_HEAD(sn, head);
    UPDATE_STORAGE(st, sh, sn->onwp[i], m, NULL);
  }
}

F_CHANGESTAT_FN(f_MultiNet){
  GET_AUX_STORAGE(StoreSubnets, sn);
  GET_STORAGE(Model*, ms);

  for(unsigned int i=1; i<=sn->ns; i++){
    if(ms[i-1]) ModelDestroy(sn->onwp[i], ms[i-1]);
  }
}

// MultiNets: Concatenate the networks' statistics; network statistic counts may be heterogeneous.

I_CHANGESTAT_FN(i_MultiNets){
  double *inputs = INPUT_PARAM; 
  GET_AUX_STORAGE(StoreSubnets, sn); inputs++;
  unsigned int ns = sn->ns;
  double *pos = inputs;
  inputs+=ns+1;
  ALLOC_STORAGE(ns, Model*, ms);

  for(unsigned int i=1; i<=sn->ns; i++){
    if(pos[i-1]!=pos[i]){
      ms[i-1] = unpack_Model_as_double(&inputs);
      InitStats(sn->onwp[i], ms[i-1]);
    }
  }
}

C_CHANGESTAT_FN(c_MultiNets){
  double *pos = INPUT_PARAM; // Starting positions of subnetworks' statistics.
  GET_AUX_STORAGE(StoreSubnets, sn); pos++;
  GET_STORAGE(Model*, ms);

  unsigned int i = MN_SID_TAIL(sn, tail);
  Vertex st = MN_IO_TAIL(sn, tail), sh = MN_IO_HEAD(sn, head);
  if(pos[i-1]!=pos[i]){
    Model *m = ms[i-1];
    ChangeStats(1, &st, &sh, sn->onwp[i], m);
    memcpy(CHANGE_STAT + (unsigned int)(pos[i-1]), m->workspace, m->n_stats*sizeof(double));
  }
}

U_CHANGESTAT_FN(u_MultiNets){
  double *pos = INPUT_PARAM; // Starting positions of subnetworks' statistics.
  GET_AUX_STORAGE(StoreSubnets, sn); pos++;
  GET_STORAGE(Model*, ms);

  unsigned int i = MN_SID_TAIL(sn, tail);
  if(pos[i-1]!=pos[i]){
    Vertex st = MN_IO_TAIL(sn, tail), sh = MN_IO_HEAD(sn, head);
    UPDATE_STORAGE(st, sh, sn->onwp[i], ms[i-1], NULL);
  }
}

F_CHANGESTAT_FN(f_MultiNets){
  double *pos = INPUT_PARAM; // Starting positions of subnetworks' statistics.
  GET_AUX_STORAGE(StoreSubnets, sn); pos++;
  GET_STORAGE(Model*, ms);

  for(unsigned int i=1; i<=sn->ns; i++){
    if(pos[i-1]!=pos[i]){
      ModelDestroy(sn->onwp[i], ms[i-1]);
    }
  }
}

// ByNetDStats

I_CHANGESTAT_FN(i_ByNetDStats){
  double *inputs = INPUT_PARAM; 
  GET_AUX_STORAGE(StoreSubnets, sn); inputs++;
  unsigned int ns = sn->ns;
  inputs+=ns+1; // Skip over subsets.

  Model *m = STORAGE = unpack_Model_as_double(&inputs);
  InitStats(nwp, m);
}

C_CHANGESTAT_FN(c_ByNetDStats){
  double *pos = INPUT_PARAM; // Starting positions of subnetworks' statistics.
  GET_AUX_STORAGE(StoreSubnets, sn); pos++;
  GET_STORAGE(Model, m);

  unsigned int i = MN_SID_TAIL(sn, tail);
  if(pos[i-1]!=pos[i]){
    ChangeStats(1, &tail, &head, nwp, m);
    memcpy(CHANGE_STAT + (unsigned int)pos[i], m->workspace, m->n_stats*sizeof(double));
  }
}

U_CHANGESTAT_FN(u_ByNetDStats){
  GET_STORAGE(Model, m);
  UPDATE_STORAGE(tail, head, nwp, m, NULL);
}

F_CHANGESTAT_FN(f_ByNetDStats){
  GET_STORAGE(Model, m);
  ModelDestroy(nwp, m);
  STORAGE = NULL;
}

