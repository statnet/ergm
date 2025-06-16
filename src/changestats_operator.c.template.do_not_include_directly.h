/*  File src/wtchangestats_operator.c in package ergm, part of the Statnet
 *  suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */

/* passthrough(formula) */

ETYPE(I_CHANGESTAT_FN)(ETYPE(i_, passthrough_term)){
  // No need to allocate it: we are only storing a pointer to a model.
  ETYPE(Model) *m = STORAGE = ETYPE(ModelInitialize)(getListElement(mtp->R, "submodel"), mtp->ext_state,  nwp, FALSE);

  ETYPE(SELECT_C_OR_D_BASED_ON_SUBMODEL)(m);
  ETYPE(DELETE_IF_UNUSED_IN_SUBMODEL)(x_func, m);
  ETYPE(DELETE_IF_UNUSED_IN_SUBMODEL)(z_func, m);
}

ETYPE(D_CHANGESTAT_FN)(ETYPE(d_, passthrough_term)){
  GET_STORAGE(ETYPE(Model), m);

  ETYPE(ChangeStats)(ntoggles, tails, heads, IFEWT(weights,) nwp, m);

  memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double));
}

ETYPE(C_CHANGESTAT_FN)(ETYPE(c_, passthrough_term)){
  GET_STORAGE(ETYPE(Model), m);

  ETYPE(ChangeStats1)(tail, head, IFEWT(weight,) nwp, m, edgestate);

  memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double));
}

ETYPE(X_CHANGESTAT_PROPAGATE_FN)(ETYPE(x_, passthrough_term), GET_STORAGE(ETYPE(Model), m), m)


ETYPE(Z_CHANGESTAT_FN)(ETYPE(z_, passthrough_term)){
  GET_STORAGE(ETYPE(Model), m);

  ETYPE(ZStats)(nwp, m, FALSE);

  memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double));
}


ETYPE(F_CHANGESTAT_FN)(ETYPE(f_, passthrough_term)){
  GET_STORAGE(ETYPE(Model), m);

  ETYPE(ModelDestroy)(nwp, m);

  STORAGE = NULL;
}


/* .submodel_and_summary(formula) */

ETYPE(I_CHANGESTAT_FN)(ETYPE(i__, submodel_and_summary_term)){
  ALLOC_AUX_STORAGE(1, ETYPE(Store, ModelAndStats), storage);

  // Unpack the submodel.
  ETYPE(Model) *m = storage->m = ETYPE(ModelInitialize)(getListElement(mtp->R, "submodel"),  NULL, nwp, FALSE);

  storage->stats = R_Calloc(m->n_stats, double);

  ETYPE(SummStats)(0, NULL, NULL, IFEWT(NULL,) nwp, m);
  memcpy(storage->stats, m->workspace, m->n_stats*sizeof(double));

  DELETE_IF_UNUSED_IN_SUBMODEL(z_func, m);
}

ETYPE(U_CHANGESTAT_FN)(ETYPE(u__, submodel_and_summary_term)){
  GET_AUX_STORAGE(ETYPE(Store, ModelAndStats), storage);
  ETYPE(Model) *m = storage->m;

  ETYPE(ChangeStats1)(tail, head, IFEWT(weight,) nwp, m, edgestate);
  addonto(storage->stats, m->workspace, m->n_stats);
}

ETYPE(F_CHANGESTAT_FN)(ETYPE(f__, submodel_and_summary_term)){
  GET_AUX_STORAGE(ETYPE(Store, ModelAndStats), storage);

  R_Free(storage->stats);
  ETYPE(ModelDestroy)(nwp, storage->m);
}


// wtSum: Take a weighted sum of the models' statistics.

ETYPE(I_CHANGESTAT_FN)(ETYPE(i_, Sum)){
  unsigned int nms = *IINPUT_PARAM;
  ALLOC_STORAGE(nms, ETYPE(Model)*, ms);

  SEXP submodels = getListElement(mtp->R, "submodels");
  for(unsigned int i=0; i<nms; i++){
    ms[i] = ETYPE(ModelInitialize)(VECTOR_ELT(submodels,i), isNULL(mtp->ext_state) ? NULL : VECTOR_ELT(mtp->ext_state,i), nwp, FALSE);
  }
  ETYPE(DELETE_IF_UNUSED_IN_SUBMODELS)(x_func, ms, nms);
  ETYPE(DELETE_IF_UNUSED_IN_SUBMODELS)(z_func, ms, nms);
}

ETYPE(C_CHANGESTAT_FN)(ETYPE(c_, Sum)){
  GET_STORAGE(ETYPE(Model)*, ms);
  unsigned int nms = *IINPUT_PARAM;
  double *wts = INPUT_PARAM;

  for(unsigned int i=0; i<nms; i++){
    ETYPE(Model) *m = ms[i];
    ETYPE(ChangeStats1)(tail, head, IFEWT(weight,) nwp, m, edgestate);
    for(unsigned int j=0; j<m->n_stats; j++)
      for(unsigned int k=0; k<N_CHANGE_STATS; k++)
	CHANGE_STAT[k] += m->workspace[j]* *(wts++);
  }
}

ETYPE(Z_CHANGESTAT_FN)(ETYPE(z_, Sum)){
  GET_STORAGE(ETYPE(Model)*, ms);
  unsigned int nms = *IINPUT_PARAM;
  double *wts = INPUT_PARAM;

  for(unsigned int i=0; i<nms; i++){
    ETYPE(Model) *m = ms[i];
    ETYPE(ZStats)(nwp, m, FALSE);
    for(unsigned int j=0; j<m->n_stats; j++)
      for(unsigned int k=0; k<N_CHANGE_STATS; k++)
	CHANGE_STAT[k] += m->workspace[j]* *(wts++);
  }
}

ETYPE(X_CHANGESTAT_FN)(ETYPE(x_, Sum)){
  GET_STORAGE(ETYPE(Model)*, ms);
  unsigned int nms = *IINPUT_PARAM;
  double *wts = INPUT_PARAM;

  for(unsigned int i=0; i<nms; i++){
    ETYPE(Model) *m = ms[i];
    ETYPE(PROPAGATE_X_SIGNAL_INTO)(nwp, m, m->workspace);
    for(unsigned int j=0; j<m->n_stats; j++)
      for(unsigned int k=0; k<N_CHANGE_STATS; k++)
	CHANGE_STAT[k] += m->workspace[j]* *(wts++);
  }
}

ETYPE(F_CHANGESTAT_FN)(ETYPE(f_, Sum)){
  GET_STORAGE(ETYPE(Model)*, ms);
  unsigned int nms = *IINPUT_PARAM;

  for(unsigned int i=0; i<nms; i++){
    ETYPE(ModelDestroy)(nwp, ms[i]);
  }
}


// Log: Take a natural logarithm of the model's statistics.

ETYPE(I_CHANGESTAT_FN)(ETYPE(i_, Log)){
  GET_AUX_STORAGE(ETYPE(Store, ModelAndStats), modstats);
  DELETE_IF_UNUSED_IN_SUBMODEL(z_func, modstats->m);
}

ETYPE(C_CHANGESTAT_FN)(ETYPE(c_, Log)){
  GET_AUX_STORAGE(ETYPE(Store, ModelAndStats), modstats);
  double *log0 = INPUT_PARAM;

  ETYPE(ChangeStats1)(tail, head, IFEWT(weight,) nwp, modstats->m, edgestate);
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

ETYPE(Z_CHANGESTAT_FN)(ETYPE(z_, Log)){
  GET_AUX_STORAGE(ETYPE(Store, ModelAndStats), modstats);
  double* log0 = INPUT_PARAM;

  ETYPE(EmptyNetworkStats)(modstats->m, FALSE);
  memcpy(CHANGE_STAT, modstats->m->workspace, N_CHANGE_STATS*sizeof(double));
  ETYPE(ZStats)(nwp, modstats->m, FALSE);
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

ETYPE(I_CHANGESTAT_FN)(ETYPE(i_, Exp)){
  GET_AUX_STORAGE(ETYPE(Store, ModelAndStats), modstats);
  DELETE_IF_UNUSED_IN_SUBMODEL(z_func, modstats->m);
}

ETYPE(C_CHANGESTAT_FN)(ETYPE(c_, Exp)){
  GET_AUX_STORAGE(ETYPE(Store, ModelAndStats), modstats);

  ETYPE(ChangeStats1)(tail, head, IFEWT(weight,) nwp, modstats->m, edgestate);
  for(unsigned int i=0; i<N_CHANGE_STATS; i++){
    if(modstats->m->workspace[i] == 0) CHANGE_STAT[i] = 0;
    else CHANGE_STAT[i] = exp(modstats->stats[i]+modstats->m->workspace[i]) - exp(modstats->stats[i]);
  }
}

ETYPE(Z_CHANGESTAT_FN)(ETYPE(z_, Exp)){
  GET_AUX_STORAGE(ETYPE(Store, ModelAndStats), modstats);
  ETYPE(EmptyNetworkStats)(modstats->m, FALSE);
  memcpy(CHANGE_STAT, modstats->m->workspace, N_CHANGE_STATS*sizeof(double));
  ETYPE(ZStats)(nwp, modstats->m, FALSE);
  for(unsigned int i=0; i<N_CHANGE_STATS; i++){
    if(modstats->m->workspace[i] == 0) CHANGE_STAT[i] = 0;
    else CHANGE_STAT[i] = exp(CHANGE_STAT[i]+modstats->m->workspace[i]) - exp(CHANGE_STAT[i]);
  }
}
