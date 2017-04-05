#include "changestat_operator.h"

/* passthrough(formula) */

I_CHANGESTAT_FN(i_passthrough_term){
  double *inputs = INPUT_PARAM;
  // No need to allocate it: we are only storing a pointer to a model.

  STORAGE = unpack_Model_as_double(&inputs);
 
  InitStats(nwp, STORAGE);
}

D_CHANGESTAT_FN(d_passthrough_term){
  GET_STORAGE(Model, m);

  ChangeStats(ntoggles, tails, heads, nwp, m);
  
  memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double));
}

U_CHANGESTAT_FN(u_passthrough_term){
  GET_STORAGE(Model, m);

  UPDATE_STORAGE(tail, head, nwp, m, NULL);
}

F_CHANGESTAT_FN(f_passthrough_term){
  GET_STORAGE(Model, m);

  ModelDestroy(nwp, m);

  STORAGE = NULL;
}

/* .submodel(formula) */

I_CHANGESTAT_FN(i__submodel_term){
  double *inputs = INPUT_PARAM+1;
  // No need to allocate it: we are only storing a pointer to a model.

  AUX_STORAGE = unpack_Model_as_double(&inputs);

  InitStats(nwp, AUX_STORAGE);
}

U_CHANGESTAT_FN(u__submodel_term){
  GET_AUX_STORAGE(Model, m);

  UPDATE_STORAGE(tail, head, nwp, m, NULL);
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
  double *inputs = INPUT_PARAM+1;
  GET_STORAGE(Model, m); // No need to allocate, since we just need a pointer.

  // Unpack the submodel.
  STORAGE = m = unpack_Model_as_double(&inputs);
  // Initialize empty network.
  Network nw = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, 0, 0, NULL);

  // Initialize storage for submodel terms 
  InitStats(&nw, m);

  ALLOC_AUX_STORAGE(m->n_stats, double, stats);
  memcpy(stats, inputs, m->n_stats*sizeof(double));

  // Evaluate the initial summary statistics (the slow way).
  for(Vertex tail=1; tail <= N_TAILS; tail++){
    Vertex head;
    Edge e;
    STEP_THROUGH_OUTEDGES(tail, e, head) {
      ChangeStats(1, &tail, &head, &nw, m);
      for(unsigned int k=0; k<m->n_stats; k++)
	stats[k] += m->workspace[k];
      UPDATE_STORAGE(tail, head, &nw, m, NULL);
      ToggleEdge(tail, head, &nw);
    }
  }
  // Note that nw is now identitical to nwp.
  NetworkDestroy(&nw);
}

U_CHANGESTAT_FN(u__summary_term){
  GET_STORAGE(Model, m);
  GET_AUX_STORAGE(double, stats);

  ChangeStats(1, &tail, &head, nwp, m);
  for(unsigned int k=0; k<m->n_stats; k++)
    stats[k] += m->workspace[k];

  UPDATE_STORAGE(tail, head, nwp, m, NULL);
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
  for(unsigned int i=0; i<INPUT_PARAM[1]; i++){
    Rprintf(" %f", stats[i]);
  }
  Rprintf(" \n");
  /* Rprintf(" ]\n"); */
}

F_CHANGESTAT_FN(f_summary_test_term){
  GET_AUX_STORAGE(double, stats);
  /* Rprintf("Test .summary auxiliary (last): ["); */
  for(unsigned int i=0; i<INPUT_PARAM[1]; i++) Rprintf(" %f", stats[i]);
  Rprintf(" \n");
  /* Rprintf(" ]\n"); */
}


