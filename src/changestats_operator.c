#include "ergm_state.h"
#include "netstats.h"
#include "ergm_changestat_operator.h"
#include "ergm_changestat_auxnet.h"
#include "ergm_util.h"

/* passthrough(formula) */

I_CHANGESTAT_FN(i_passthrough_term){
  // No need to allocate it: we are only storing a pointer to a model.
  Model *m = STORAGE = ModelInitialize(getListElement(mtp->R, "submodel"), mtp->ext_state,  nwp, FALSE);

  SELECT_C_OR_D_BASED_ON_SUBMODEL(m);
  DELETE_IF_UNUSED_IN_SUBMODEL(u_func, m);
}

D_CHANGESTAT_FN(d_passthrough_term){
  GET_STORAGE(Model, m);

  ChangeStats(ntoggles, tails, heads, nwp, m);
  
  memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double));
}

C_CHANGESTAT_FN(c_passthrough_term){
  GET_STORAGE(Model, m);

  ChangeStats1(tail, head, nwp, m, edgeflag);

  memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double));
}


U_CHANGESTAT_FN(u_passthrough_term){
  GET_STORAGE(Model, m);

  UPDATE_STORAGE(tail, head, nwp, m, NULL, edgeflag);
}

F_CHANGESTAT_FN(f_passthrough_term){
  GET_STORAGE(Model, m);

  ModelDestroy(nwp, m);

  STORAGE = NULL;
}

/* .submodel(formula) */

I_CHANGESTAT_FN(i__submodel_term){
  // No need to allocate it: we are only storing a pointer to a model.
  Model *m = AUX_STORAGE = ModelInitialize(getListElement(mtp->R, "submodel"), mtp->ext_state,  nwp, FALSE);
  DELETE_IF_UNUSED_IN_SUBMODEL(u_func, m);
}

U_CHANGESTAT_FN(u__submodel_term){
  GET_AUX_STORAGE(Model, m);

  UPDATE_STORAGE(tail, head, nwp, m, NULL, edgeflag);
}

F_CHANGESTAT_FN(f__submodel_term){
  GET_AUX_STORAGE(Model, m);

  ModelDestroy(nwp, m);

  AUX_STORAGE=NULL;
}

/* submodel.test(formula) */

D_CHANGESTAT_FN(d_submodel_test_term){
  GET_AUX_STORAGE(Model, m);

  double *tmp = m->workspace;
  m->workspace = CHANGE_STAT;
  ChangeStats(ntoggles, tails, heads, nwp, m);
  m->workspace = tmp;
}

/* .summary(formula) */


I_CHANGESTAT_FN(i__summary_term){
  GET_STORAGE(Model, m); // No need to allocate, since we just need a pointer.

  // Initialize empty network.
  Network *tmpnwp = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, 0, 0, NULL);
  // Unpack the submodel.
  STORAGE = m = ModelInitialize(getListElement(mtp->R, "submodel"), mtp->ext_state,  tmpnwp, FALSE);

  ALLOC_AUX_STORAGE(m->n_stats, double, stats);
  ErgmState s={.stats=NULL,
               .nwp=tmpnwp,
               .m=m,
               .MHp=NULL};
  Vertex *tails = Calloc(EDGECOUNT(nwp), Vertex);
  Vertex *heads = Calloc(EDGECOUNT(nwp), Vertex);
  EdgeTree2EdgeList(tails, heads, nwp, EDGECOUNT(nwp));

  SummStats(&s, EDGECOUNT(nwp), tails, heads, stats);
  // Note that nw is now identitical to nwp.
  NetworkDestroy(tmpnwp);
}

U_CHANGESTAT_FN(u__summary_term){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(double, stats);

  ChangeStats(1, &tail, &head, nwp, m);
  addonto(stats, m->workspace, m->n_stats);

  UPDATE_STORAGE(tail, head, nwp, m, NULL, edgeflag);
}

F_CHANGESTAT_FN(f__summary_term){
  GET_STORAGE(Model, m);

  ModelDestroy(nwp, m);

  STORAGE=NULL;
}

/* summary.test(formula) */

C_CHANGESTAT_FN(c_summary_test_term){
  GET_AUX_STORAGE(double, stats);
  CHANGE_STAT[0] = 0;
  /* Rprintf("Test .summary auxiliary: ["); */
  for(unsigned int i=0; i<INPUT_PARAM[0]; i++){
    Rprintf(" %f", stats[i]);
  }
  Rprintf(" \n");
  /* Rprintf(" ]\n"); */
}

F_CHANGESTAT_FN(f_summary_test_term){
  GET_AUX_STORAGE(double, stats);
  /* Rprintf("Test .summary auxiliary (last): ["); */
  for(unsigned int i=0; i<INPUT_PARAM[0]; i++) Rprintf(" %f", stats[i]);
  Rprintf(" \n");
  /* Rprintf(" ]\n"); */
}

// Sum: Take a weighted sum of the models' statistics.

I_CHANGESTAT_FN(i_Sum){
  /*
    inputs expected:
    1: number of models (nms)
    1: total length of all weight matrices (tml)
    tml: a list of mapping matrices in row-major order
  */
  
  double *inputs = INPUT_PARAM; 
  unsigned int nms = *(inputs++);
  unsigned int tml = *(inputs++);
  inputs+=tml; // Jump to start of model specifications.
  
  ALLOC_STORAGE(nms, Model*, ms);

  SEXP submodels = getListElement(mtp->R, "submodels");
  for(unsigned int i=0; i<nms; i++){
    ms[i] = ModelInitialize(VECTOR_ELT(submodels, i), isNULL(mtp->ext_state) ? NULL : VECTOR_ELT(mtp->ext_state,i), nwp, FALSE);
  }
  DELETE_IF_UNUSED_IN_SUBMODELS(u_func, ms, nms);
}

C_CHANGESTAT_FN(c_Sum){
  double *inputs = INPUT_PARAM; 
  GET_STORAGE(Model*, ms);
  unsigned int nms = *(inputs++);
  inputs++; //  Skip total length of weight matrices.
  double *wts = inputs;

  for(unsigned int i=0; i<nms; i++){
    Model *m = ms[i];
    ChangeStats(1, &tail, &head, nwp, m);
    for(unsigned int j=0; j<m->n_stats; j++)
      for(unsigned int k=0; k<N_CHANGE_STATS; k++)
	CHANGE_STAT[k] += m->workspace[j]* *(wts++);
  }
}

U_CHANGESTAT_FN(u_Sum){
  double *inputs = INPUT_PARAM; 
  GET_STORAGE(Model*, ms);
  unsigned int nms = *(inputs++);

  for(unsigned int i=0; i<nms; i++){
    Model *m = ms[i];
    UPDATE_STORAGE(tail, head, nwp, m, NULL, edgeflag);
  }
}

F_CHANGESTAT_FN(f_Sum){
  double *inputs = INPUT_PARAM; 
  GET_STORAGE(Model*, ms);
  unsigned int nms = *(inputs++);

  for(unsigned int i=0; i<nms; i++){
    ModelDestroy(nwp, ms[i]);
  }
}

#include "ergm_changestats_auxnet.h"

ON_AUXNET(_discord_net_Network)
ON_AUXNET(_intersect_net_Network)
ON_AUXNET(_union_net_Network)
ON_AUXNET(_blockdiag_net)
ON_AUXNET(_undir_net)
ON_AUXNET(_filter_formula_net)
ON_AUXNET(_subgraph_net)
