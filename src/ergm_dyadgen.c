#include "ergm_dyadgen.h"
#include "ergm_Rutil.h"

DyadGen *DyadGenInitialize(enum DyadGenType type, void *data){
  DyadGen *gen = Calloc(1, DyadGen);
  gen->type = type;

  switch(gen->type){
  case RandDyadGen:
    gen->data.nwp = data;
    gen->ndyads = DYADCOUNT(gen->data.nwp);
    break;
  case WtRandDyadGen:
    gen->data.wtnwp = data;
    gen->ndyads = DYADCOUNT(gen->data.wtnwp);
    break;
  case RLEBDM1DGen:
    // data must be a double ** that is shifted.
    gen->data.rlebdm = unpack_RLEBDM1D(data);
    gen->ndyads = gen->data.rlebdm.ndyads;
    break;
  case EdgeListGen:
    // data must be an int ** that is shifted.
    gen->data.el = *(int **)data;
    gen->ndyads = *gen->data.el;
    (*(int **)data) += gen->ndyads*2 + 1;
    break;
  default:
    error("Undefined dyad generator type.");
  }

  return gen;
}


DyadGen *DyadGenInitializeR(SEXP pR, void *any_nwp){
  // If there is a dyadgen element, use that; otherwise assume we're already getting a dyadgen list.
  SEXP dgR = getListElement(pR, "dyadgen");
  if(isNULL(dgR)) dgR = pR;

  enum DyadGenType type = asInteger(getListElement(dgR, "type"));

  DyadGen *gen = Calloc(1, DyadGen);
  gen->type = type;

  switch(gen->type){
  case RandDyadGen:
    return DyadGenInitialize(type, any_nwp);
    break;
  case WtRandDyadGen:
    return DyadGenInitialize(type, any_nwp);
    break;
  case RLEBDM1DGen:
    // RLEBDM1D's unpacking function expects a double **.
    {
      double *tmp = REAL(getListElement(pR, "data"));
      return DyadGenInitialize(type, &tmp);
    }
    break;
  case EdgeListGen:
    // RLEBDM1D's unpacking function expects an int **.
    {
      int *tmp = INTEGER(getListElement(dgR, "data"));
      return DyadGenInitialize(type, &tmp);
    }
    break;
  default:
    error("Undefined dyad generator type.");
  }
}


void DyadGenDestroy(DyadGen *gen){
  Free(gen);
}
