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

#include "ergm_edgetype_set_binary.h"
#include "changestats_operator.c.template.do_not_include_directly.h"

#include "ergm_changestats_auxnet.h"

ON_AUXNET(_discord_net_Network)
ON_AUXNET(_intersect_net_Network)
ON_AUXNET(_union_net_Network)
ON_AUXNET(_blockdiag_net)
ON_AUXNET(_undir_net)
ON_AUXNET(_filter_formula_net)
ON_AUXNET(_subgraph_net)
