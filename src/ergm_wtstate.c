#include "ergm_wtstate.h"
#include "ergm_constants.h"

static WtErgmState **ergm_wtstate_array = NULL;
static unsigned int ergm_wtstate_array_len = 0;
static unsigned int ergm_wtstate_array_maxlen = 0;

WtErgmState *WtErgmStateInit(SEXP stateR,
                             unsigned int flags){
  WtErgmState *s = Calloc(1, WtErgmState);

  /* Extract stats vector */
  SEXP tmp = getListElement(stateR, "stats");
  s->stats = length(tmp) ? REAL(tmp) : NULL;

  /* Form the network */
  s->nwp=Redgelist2WtNetwork(getListElement(stateR,"el"), flags & ERGM_STATE_EMPTY_NET);

  /* Initialize the model */
  s->m=NULL;
  tmp = getListElement(stateR, "model");
  if(s->nwp && length(tmp)){ // Model also requires network.
    if(asInteger(getListElement(stateR, "ext.flag"))==ERGM_STATE_R_CHANGED) error("R ergm_state has changed in R but has not been reconciled.");
    s->m = WtModelInitialize(tmp, getListElement(stateR, "ext.state"), s->nwp, flags & ERGM_STATE_NO_INIT_S);
  }

  /* Initialize the M-H proposal */
  s->MHp=NULL;
  if(!(flags & ERGM_STATE_NO_INIT_PROP) && s->m && length(tmp = getListElement(stateR, "proposal"))) // Proposal also requires model's auxiliaries.
    s->MHp = WtMHProposalInitialize(tmp, s->nwp, s->m->termarray->aux_storage);

  if(ergm_wtstate_array_len == ergm_wtstate_array_maxlen){
    ergm_wtstate_array_maxlen = MAX(1, ergm_wtstate_array_maxlen*2);
    ergm_wtstate_array = Realloc(ergm_wtstate_array, ergm_wtstate_array_maxlen, WtErgmState*);
  }
  ergm_wtstate_array[ergm_wtstate_array_len++] = s;

  return s;
}

SEXP WtErgmStateRSave(SEXP startR, WtErgmState *s){
  // Duplicate state
  SEXP outl = PROTECT(allocVector(VECSXP, length(startR)));
  setAttrib(outl, R_NamesSymbol, getAttrib(startR, R_NamesSymbol));
  for(unsigned int i=0; i<length(startR); i++)
    SET_VECTOR_ELT(outl, i, VECTOR_ELT(startR, i));

  // Network state
  if(s->nwp) setListElement(outl, "el", WtNetwork2Redgelist(s->nwp));

  // Extended state
  if(s->m){{ // To limit the scope of the variables.
      SEXP ext_l = PROTECT(allocVector(VECSXP, s->m->n_terms));
      unsigned int i=0;
      WtEXEC_THROUGH_TERMS(s->m, {
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

void WtErgmStateDestroy(WtErgmState *s){
  // Find and clear the corresponding element in the active state
  // array. Note that most of the time, this will be the first element
  // in the array.
  unsigned int i=0;
  while(ergm_wtstate_array[i] != s) i++;
  ergm_wtstate_array[i] = ergm_wtstate_array[--ergm_wtstate_array_len];
  ergm_wtstate_array[ergm_wtstate_array_len] = NULL;

  if(s->MHp) WtMHProposalDestroy(s->MHp, s->nwp);
  if(s->m) WtModelDestroy(s->nwp, s->m);
  if(s->nwp) WtNetworkDestroy(s->nwp);
  Free(s);
}

SEXP ErgmWtStateArrayClear(){
  while(ergm_wtstate_array_len) WtErgmStateDestroy(ergm_wtstate_array[0]);
  ergm_wtstate_array_maxlen = 0;
  Free(ergm_wtstate_array);
  return R_NilValue;
}
