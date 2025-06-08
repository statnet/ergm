/*  File src/changestats_operator.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#include "ergm_changestat_operator.h"
#include "ergm_changestat_auxnet.h"
#include "ergm_changestats_operator.h"
#include "ergm_util.h"

/* passthrough(formula) */

I_CHANGESTAT_FN(i_passthrough_term){
  // No need to allocate it: we are only storing a pointer to a model.
  Model *m = STORAGE = ModelInitialize(getListElement(mtp->R, "submodel"), mtp->ext_state,  nwp, FALSE);

  SELECT_C_OR_D_BASED_ON_SUBMODEL(m);
  DELETE_IF_UNUSED_IN_SUBMODEL(x_func, m);
  DELETE_IF_UNUSED_IN_SUBMODEL(z_func, m);
}

D_CHANGESTAT_FN(d_passthrough_term){
  GET_STORAGE(Model, m);

  ChangeStats(ntoggles, tails, heads, nwp, m);
  
  memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double));
}

C_CHANGESTAT_FN(c_passthrough_term){
  GET_STORAGE(Model, m);

  ChangeStats1(tail, head, nwp, m, edgestate);

  memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double));
}

X_CHANGESTAT_PROPAGATE_FN(x_passthrough_term, GET_STORAGE(Model, m), m)

Z_CHANGESTAT_FN(z_passthrough_term){
  GET_STORAGE(Model, m);

  ZStats(nwp, m, FALSE);
  memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double));
}

F_CHANGESTAT_FN(f_passthrough_term){
  GET_STORAGE(Model, m);

  ModelDestroy(nwp, m);

  STORAGE = NULL;
}

/* .submodel(formula) */

I_CHANGESTAT_FN(i__submodel_term){
  // No need to allocate it: we are only storing a pointer to a model.
  AUX_STORAGE = ModelInitialize(getListElement(mtp->R, "submodel"), NULL,  nwp, FALSE);
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
  // Unpack the submodel.
  STORAGE = m = ModelInitialize(getListElement(mtp->R, "submodel"), mtp->ext_state,  nwp, FALSE);
  ALLOC_AUX_STORAGE(m->n_stats, double, stats);

  SummStats(0, NULL, NULL, nwp, m);
  memcpy(stats, m->workspace, m->n_stats*sizeof(double));

  DELETE_IF_UNUSED_IN_SUBMODEL(z_func, m);
}

X_CHANGESTAT_FN(x__summary_term){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(double, stats);

  PROPAGATE_X_SIGNAL_ADDONTO(nwp, m, stats);
}

U_CHANGESTAT_FN(u__summary_term){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(double, stats);

  ChangeStats1(tail, head, nwp, m, edgestate);
  addonto(stats, m->workspace, m->n_stats);
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


/* .submodel_and_summary(formula) */

I_CHANGESTAT_FN(i__submodel_and_summary_term){
  ALLOC_AUX_STORAGE(1, StoreModelAndStats, storage);

  // Unpack the submodel.
  Model *m = storage->m = ModelInitialize(getListElement(mtp->R, "submodel"),  NULL, nwp, FALSE);

  storage->stats = R_Calloc(m->n_stats, double);

  SummStats(0, NULL, NULL, nwp, m);
  memcpy(storage->stats, m->workspace, m->n_stats*sizeof(double));

  DELETE_IF_UNUSED_IN_SUBMODEL(z_func, m);
}

U_CHANGESTAT_FN(u__submodel_and_summary_term){
  GET_AUX_STORAGE(StoreModelAndStats, storage);
  Model *m = storage->m;

  ChangeStats1(tail, head, nwp, m, edgestate);
  addonto(storage->stats, m->workspace, m->n_stats);
}

F_CHANGESTAT_FN(f__submodel_and_summary_term){
  GET_AUX_STORAGE(StoreModelAndStats, storage);

  R_Free(storage->stats);
  ModelDestroy(nwp, storage->m);
}


// Sum: Take a weighted sum of the models' statistics.

I_CHANGESTAT_FN(i_Sum){
  unsigned int nms = *IINPUT_PARAM;
  ALLOC_STORAGE(nms, Model*, ms);

  SEXP submodels = getListElement(mtp->R, "submodels");
  for(unsigned int i=0; i<nms; i++){
    ms[i] = ModelInitialize(VECTOR_ELT(submodels, i), isNULL(mtp->ext_state) ? NULL : VECTOR_ELT(mtp->ext_state,i), nwp, FALSE);
  }
  DELETE_IF_UNUSED_IN_SUBMODELS(x_func, ms, nms);
  DELETE_IF_UNUSED_IN_SUBMODELS(z_func, ms, nms);
}

C_CHANGESTAT_FN(c_Sum){
  GET_STORAGE(Model*, ms);
  unsigned int nms = *IINPUT_PARAM;
  double *wts = INPUT_PARAM;

  for(unsigned int i=0; i<nms; i++){
    Model *m = ms[i];
    ChangeStats1(tail, head, nwp, m, edgestate);
    for(unsigned int j=0; j<m->n_stats; j++)
      for(unsigned int k=0; k<N_CHANGE_STATS; k++)
	CHANGE_STAT[k] += m->workspace[j]* *(wts++);
  }
}

Z_CHANGESTAT_FN(z_Sum){
  GET_STORAGE(Model*, ms);
  unsigned int nms = *IINPUT_PARAM;
  double *wts = INPUT_PARAM;

  for(unsigned int i=0; i<nms; i++){
    Model *m = ms[i];
    ZStats(nwp, m, FALSE);
    for(unsigned int j=0; j<m->n_stats; j++)
      for(unsigned int k=0; k<N_CHANGE_STATS; k++)
	CHANGE_STAT[k] += m->workspace[j]* *(wts++);
  }
}

X_CHANGESTAT_FN(x_Sum){
  GET_STORAGE(Model*, ms);
  unsigned int nms = *IINPUT_PARAM;
  double *wts = INPUT_PARAM;

  for(unsigned int i=0; i<nms; i++){
    Model *m = ms[i];
    PROPAGATE_X_SIGNAL_INTO(nwp, m, m->workspace);
    for(unsigned int j=0; j<m->n_stats; j++)
      for(unsigned int k=0; k<N_CHANGE_STATS; k++)
	CHANGE_STAT[k] += m->workspace[j]* *(wts++);
  }
}

F_CHANGESTAT_FN(f_Sum){
  GET_STORAGE(Model*, ms);
  unsigned int nms = *IINPUT_PARAM;

  for(unsigned int i=0; i<nms; i++){
    ModelDestroy(nwp, ms[i]);
  }
}


// Log: Take a natural logarithm of the model's statistics.

I_CHANGESTAT_FN(i_Log){
  GET_AUX_STORAGE(StoreModelAndStats, modstats);
  DELETE_IF_UNUSED_IN_SUBMODEL(z_func, modstats->m);
}

C_CHANGESTAT_FN(c_Log){
  GET_AUX_STORAGE(StoreModelAndStats, modstats);
  double *log0 = INPUT_PARAM;

  ChangeStats1(tail, head, nwp, modstats->m, edgestate);
  for(unsigned int i=0; i<N_CHANGE_STATS; i++){
    if(modstats->m->workspace[i] == 0) CHANGE_STAT[i] = 0;
    else{
      double old = modstats->stats[i];
      old = old == 0 ? log0[i] : log(old);
      double new = modstats->stats[i]+modstats->m->workspace[i];
      new = new == 0 ? log0[i] : log(new);
      CHANGE_STAT[i] = new - old;
    }
  }
}

Z_CHANGESTAT_FN(z_Log){
  GET_AUX_STORAGE(StoreModelAndStats, modstats);
  double* log0 = INPUT_PARAM;

  EmptyNetworkStats(modstats->m, FALSE);
  memcpy(CHANGE_STAT, modstats->m->workspace, N_CHANGE_STATS*sizeof(double));
  ZStats(nwp, modstats->m, FALSE);
  for(unsigned int i=0; i<N_CHANGE_STATS; i++){
    if(modstats->m->workspace[i] == 0) CHANGE_STAT[i] = 0;
    else{
      double old = CHANGE_STAT[i];
      old = old == 0 ? log0[i] : log(old);
      double new = CHANGE_STAT[i]+modstats->m->workspace[i];
      new = new == 0 ? log0[i] : log(new);
      CHANGE_STAT[i] = new - old;
    }
  }
}

// Exp: Exponentiate the model's statistics.

I_CHANGESTAT_FN(i_Exp){
  GET_AUX_STORAGE(StoreModelAndStats, modstats);
  DELETE_IF_UNUSED_IN_SUBMODEL(z_func, modstats->m);
}

C_CHANGESTAT_FN(c_Exp){
  GET_AUX_STORAGE(StoreModelAndStats, modstats);

  ChangeStats1(tail, head, nwp, modstats->m, edgestate);
  for(unsigned int i=0; i<N_CHANGE_STATS; i++){
    if(modstats->m->workspace[i] == 0) CHANGE_STAT[i] = 0;
    else CHANGE_STAT[i] = exp(modstats->stats[i]+modstats->m->workspace[i]) - exp(modstats->stats[i]);
  }
}

Z_CHANGESTAT_FN(z_Exp){
  GET_AUX_STORAGE(StoreModelAndStats, modstats);
  EmptyNetworkStats(modstats->m, FALSE);
  memcpy(CHANGE_STAT, modstats->m->workspace, N_CHANGE_STATS*sizeof(double));
  ZStats(nwp, modstats->m, FALSE);
  for(unsigned int i=0; i<N_CHANGE_STATS; i++){
    if(modstats->m->workspace[i] == 0) CHANGE_STAT[i] = 0;
    else CHANGE_STAT[i] = exp(CHANGE_STAT[i]+modstats->m->workspace[i]) - exp(CHANGE_STAT[i]);
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
