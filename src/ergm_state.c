#include "ergm_state.h"
#include "ergm_constants.h"

ErgmState *ErgmStateInit(SEXP stateR,
                         Rboolean empty, Rboolean noinit_s){
  ErgmState *s = Calloc(1, ErgmState);

  /* Extract stats vector */
  SEXP tmp = getListElement(stateR, "stats");
  s->stats = length(tmp) ? REAL(tmp) : NULL;

  /* Form the network */
  s->nwp=Redgelist2Network(getListElement(stateR,"el"), empty);

  /* Initialize the model */
  s->m=NULL;
  tmp = getListElement(stateR, "model");
  if(s->nwp && length(tmp)){ // Model also requires network.
    if(asInteger(getListElement(stateR, "ext.flag"))==ERGM_STATE_R_CHANGED) error("R ergm_state has changed in R but has not been reconciled.");
    s->m = ModelInitialize(tmp, getListElement(stateR, "ext.state"), s->nwp, noinit_s);
  }

  /* Initialize the M-H proposal */
  s->MHp=NULL;
  tmp = getListElement(stateR, "proposal");
  if(s->m && length(tmp)) // Proposal also requires model's auxiliaries.
    s->MHp = MHProposalInitialize(tmp, s->nwp, s->m->termarray->aux_storage);

  return s;
}

SEXP ErgmStateRSave(SEXP startR, ErgmState *s){
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
  if(s->MHp) MHProposalDestroy(s->MHp, s->nwp);
  if(s->m) ModelDestroy(s->nwp, s->m);
  if(s->nwp) NetworkDestroy(s->nwp);
  Free(s);
}
