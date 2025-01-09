/*  File src/ergm_state.c in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#include "ergm_kvec.h"
#include "ergm_state.h"
#include "ergm_constants.h"

static kvec_t(ErgmState *) ergm_state_array = kv_blank;

ErgmState *ErgmStateInit(SEXP stateR,
                         unsigned int flags){
  ErgmState *s = R_Calloc(1, ErgmState);

  /* Save a reference to the corresponding R object */
  s->R = stateR;

  /* Extract stats vector */
  SEXP tmp = getListElement(stateR, "stats");
  s->stats = length(tmp) ? REAL(tmp) : NULL;

  /* Form the network */
  s->nwp=Redgelist2Network(getListElement(stateR,"el"), flags & ERGM_STATE_EMPTY_NET);

  /* Initialize the model */
  s->m=NULL;
  tmp = getListElement(stateR, "model");
  if(s->nwp && length(tmp)){ // Model also requires network.
    if(asInteger(getListElement(stateR, "ext.flag"))==ERGM_STATE_R_CHANGED) error("R ergm_state has changed in R but has not been reconciled.");
    s->m = ModelInitialize(tmp, getListElement(stateR, "ext.state"), s->nwp, flags & ERGM_STATE_NO_INIT_S);
  }

  /* Initialize the M-H proposal */
  s->MHp=NULL;
  if(!(flags & ERGM_STATE_NO_INIT_PROP) && s->m && length(tmp = getListElement(stateR, "proposal"))) // Proposal also requires model's auxiliaries.
    s->MHp = MHProposalInitialize(tmp, s->nwp, s->m->termarray->aux_storage);

  kv_push(ErgmState*, ergm_state_array, s);

  return s;
}

SEXP ErgmStateRSave(ErgmState *s){
  SEXP startR = s->R;

  // Duplicate state
  SEXP outl = PROTECT(allocVector(VECSXP, length(startR)));
  setAttrib(outl, R_NamesSymbol, getAttrib(startR, R_NamesSymbol));
  for(unsigned int i=0; i<length(startR); i++)
    SET_VECTOR_ELT(outl, i, VECTOR_ELT(startR, i));

  // Network state
  if(s->nwp) setListElement(outl, "el", Network2Redgelist(s->nwp));

  // Extended state
  if(s->m){{ // To limit the scope of the variables.
      SEXP ext_l = PROTECT(allocVector(VECSXP, s->m->n_terms));
      unsigned int i=0;
      EXEC_THROUGH_TERMS(s->m, {
          if(mtp->w_func) SET_VECTOR_ELT(ext_l, i, mtp->w_func(mtp, s->nwp));
          i++;
        });
      setListElement(outl, "ext.state", ext_l);
      setListElement(outl, "ext.flag", ScalarInteger(ERGM_STATE_C_CHANGED));
      UNPROTECT(1);
    }}
  
  // Statistics
  if(s->stats){
    SEXP statsR = PROTECT(allocVector(REALSXP, length(getListElement(startR, "stats"))));
    memcpy(REAL(statsR), s->stats, length(statsR)*sizeof(double));
    setListElement(outl, "stats", statsR);
    UNPROTECT(1); // statsR
  }

  classgets(outl, getAttrib(startR, R_ClassSymbol));

  UNPROTECT(1); // outl
  return outl;
}

void ErgmStateDestroy(ErgmState *s){
  // Find and clear the corresponding element in the active state
  // array. Note that most of the time, this will be the first element
  // in the array.
  unsigned int i=0;
  while(kv_A(ergm_state_array, i) != s) i++;
  kv_del_plug(ergm_state_array, i);

  if(s->MHp) MHProposalDestroy(s->MHp, s->nwp);
  if(s->m) ModelDestroy(s->nwp, s->m);
  if(s->nwp) NetworkDestroy(s->nwp);
  R_Free(s);
}

SEXP ErgmStateArrayClear(void){
  while(kv_size(ergm_state_array)) ErgmStateDestroy(kv_A(ergm_state_array, 0));
  kv_destroy(ergm_state_array);
  return R_NilValue;
}
