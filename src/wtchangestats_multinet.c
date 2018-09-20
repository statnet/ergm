#include "ergm_wtchangestat_multinet.h"
#include "ergm_wtchangestat.h"
#include "ergm_wtchangestat_operator.h"


WtI_CHANGESTAT_FN(i__wtsubnets){
  double *inputs = INPUT_PARAM;
  ALLOC_AUX_STORAGE(1, StoreWtSubnets, sn); inputs++;
  sn->ns = *(inputs++);
  sn->inwp = nwp;
  sn->onwp = Calloc(sn->ns, WtNetwork *);
  sn->onwp--; // The -- is because WtNetwork IDs count from 1.

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

    sn->onwp[i] = WtNetworkInitialize(NULL, NULL, NULL, 0, lnnodes, DIRECTED, lbip, 0, 0, NULL);
  }
  
  EXEC_THROUGH_NET_EDGES_PRE(t, h, e, w, {
      WtSetEdge(MN_IO_TAIL(sn, t), MN_IO_HEAD(sn, h), w, sn->onwp[MN_SID_TAIL(sn, t)]);
    });
}

WtU_CHANGESTAT_FN(u__wtsubnets){ 
  GET_AUX_STORAGE(StoreWtSubnets, sn);
  WtSetEdge(MN_IO_TAIL(sn, tail), MN_IO_HEAD(sn, head), weight, sn->onwp[MN_SID_TAIL(sn, tail)]);
}

WtF_CHANGESTAT_FN(f__wtsubnets){
  GET_AUX_STORAGE(StoreWtSubnets, sn);
  for(unsigned int i=1; i<=sn->ns; i++)
    WtNetworkDestroy(sn->onwp[i]);
  sn->onwp++;
  Free(sn->onwp);
}

// MultiNet: Take a weighted networkwise sum of the networks' statistics.

WtI_CHANGESTAT_FN(i_wtMultiNet){
  /*
    inputs expects:
    1: position of the subnets auxiliary
    1: number of weights per network (nwts)
    nwts*ns: matrix of weights, in network-major order
    ?*ns: submodel specifications for nm submodels
  */
  
  double *inputs = INPUT_PARAM; 
  GET_AUX_STORAGE(StoreWtSubnets, sn); inputs++;
  unsigned int ns = sn->ns;
  unsigned int nwts = *(inputs++);
  double *wts = inputs; inputs+=ns*nwts;
  
  ALLOC_STORAGE(ns, WtModel*, ms);

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
      ms[i-1] = unpack_WtModel_as_double(&inputs);
      WtInitStats(sn->onwp[i], ms[i-1]);
    }else ms[i-1] = NULL;
  }
}

WtC_CHANGESTAT_FN(c_wtMultiNet){
  double *inputs = INPUT_PARAM; 
  GET_AUX_STORAGE(StoreWtSubnets, sn); inputs++;
  GET_STORAGE(WtModel*, ms);
  unsigned int nwts = *(inputs++);
  double *wts = inputs;

  unsigned int i = MN_SID_TAIL(sn, tail);
  WtModel *m = ms[i-1];
  if(m){ // NULL if network has weights 0.
    Vertex st = MN_IO_TAIL(sn, tail), sh = MN_IO_HEAD(sn, head);
    WtChangeStats(1, &st, &sh, &weight, sn->onwp[i], m);

    wts += (i-1)*nwts; // Position of that network's weight vector.
    for(unsigned int j=0; j<m->n_stats; j++)
      for(unsigned int k=0; k<nwts; k++)
	CHANGE_STAT[j*nwts+k] += m->workspace[j]*wts[k];
  }
}

WtU_CHANGESTAT_FN(u_wtMultiNet){
  GET_AUX_STORAGE(StoreWtSubnets, sn);
  GET_STORAGE(WtModel*, ms);

  unsigned int i = MN_SID_TAIL(sn, tail);
  WtModel *m = ms[i-1];
  if(m){ // NULL if network has weights 0.
    Vertex st = MN_IO_TAIL(sn, tail), sh = MN_IO_HEAD(sn, head);
    WtUPDATE_STORAGE(st, sh, weight, sn->onwp[i], m, NULL);
  }
}

WtF_CHANGESTAT_FN(f_wtMultiNet){
  GET_AUX_STORAGE(StoreWtSubnets, sn);
  GET_STORAGE(WtModel*, ms);

  for(unsigned int i=1; i<=sn->ns; i++){
    if(ms[i-1]) WtModelDestroy(sn->onwp[i], ms[i-1]);
  }
}

// MultiNets: Concatenate the networks' statistics; network statistic counts may be heterogeneous.

WtI_CHANGESTAT_FN(i_wtMultiNets){
  double *inputs = INPUT_PARAM; 
  GET_AUX_STORAGE(StoreWtSubnets, sn); inputs++;
  unsigned int ns = sn->ns;
  double *pos = inputs;
  inputs+=ns+1;
  ALLOC_STORAGE(ns, WtModel*, ms);

  for(unsigned int i=1; i<=sn->ns; i++){
    if(pos[i-1]!=pos[i]){
      ms[i-1] = unpack_WtModel_as_double(&inputs);
      WtInitStats(sn->onwp[i], ms[i-1]);
    }
  }
}

WtC_CHANGESTAT_FN(c_wtMultiNets){
  double *pos = INPUT_PARAM; // Starting positions of subnetworks' statistics.
  GET_AUX_STORAGE(StoreWtSubnets, sn); pos++;
  GET_STORAGE(WtModel*, ms);

  unsigned int i = MN_SID_TAIL(sn, tail);
  Vertex st = MN_IO_TAIL(sn, tail), sh = MN_IO_HEAD(sn, head);
  if(pos[i-1]!=pos[i]){
    WtModel *m = ms[i-1];
    WtChangeStats(1, &st, &sh, &weight, sn->onwp[i], m);
    memcpy(CHANGE_STAT + (unsigned int)(pos[i-1]), m->workspace, m->n_stats*sizeof(double));
  }
}

WtU_CHANGESTAT_FN(u_wtMultiNets){
  double *pos = INPUT_PARAM; // Starting positions of subnetworks' statistics.
  GET_AUX_STORAGE(StoreWtSubnets, sn); pos++;
  GET_STORAGE(WtModel*, ms);

  unsigned int i = MN_SID_TAIL(sn, tail);
  if(pos[i-1]!=pos[i]){
    Vertex st = MN_IO_TAIL(sn, tail), sh = MN_IO_HEAD(sn, head);
    WtUPDATE_STORAGE(st, sh, weight, sn->onwp[i], ms[i-1], NULL);
  }
}

WtF_CHANGESTAT_FN(f_wtMultiNets){
  double *pos = INPUT_PARAM; // Starting positions of subnetworks' statistics.
  GET_AUX_STORAGE(StoreWtSubnets, sn); pos++;
  GET_STORAGE(WtModel*, ms);

  for(unsigned int i=1; i<=sn->ns; i++){
    if(pos[i-1]!=pos[i]){
      WtModelDestroy(sn->onwp[i], ms[i-1]);
    }
  }
}

