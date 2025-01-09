/*  File src/ergm_dyadgen.c in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#include "ergm_dyadgen.h"
#include "ergm_Rutil.h"
#include "ergm_changestat.h"
#include "ergm_wtchangestat.h"
#include "ergm_kvec.h"

typedef struct {
  OnDyadGenInit callback;
  void *payload;
} OnDyadGenInitInfo;

static kvec_t(OnDyadGenInitInfo) dyadgen_init_callbacks = kv_blank;

void DyadGenSetUpIntersect(DyadGen *gen, void *track_nwp, Rboolean force){
    switch(gen->type){
    case RandDyadGen:
    case WtRandDyadGen:
      gen->inter.uel = NULL;
      break;
    case EdgeListGen:
    case RLEBDM1DGen:
      {
        Network *nwp = track_nwp ? track_nwp : gen->nwp.b;
        gen->nwp.b = nwp;
        gen->inter.uel = UnsrtELInitialize(0, NULL, NULL, FALSE);
        EXEC_THROUGH_NET_EDGES(t, h, e, {
            if(DyadGenSearch(t, h, gen)){
              UnsrtELInsert(t, h, gen->inter.uel);
            }
          });
        
        if(!force && UnsrtELSize(gen->inter.uel) == EDGECOUNT(nwp)){ // There are no ties in the initial network that are fixed.
          UnsrtELDestroy(gen->inter.uel);
          gen->inter.uel = NULL; // "Signal" that there is no discordance network.
        }else{
          gen->intertype = UnsrtELDyadGen;
          AddOnNetworkEdgeChange(nwp, (OnNetworkEdgeChange) DyadGenUpdate, gen, INT_MAX);
        }
      }
      break;
    case WtEdgeListGen:
    case WtRLEBDM1DGen:
      {
        WtNetwork *nwp = track_nwp ? track_nwp : gen->nwp.w;
        gen->nwp.w = nwp;
        gen->inter.uel = UnsrtELInitialize(0, NULL, NULL, FALSE);
        WtEXEC_THROUGH_NET_EDGES(t, h, w, e, {
            (void) e;
            if(w!=0 && DyadGenSearch(t, h, gen)){
              UnsrtELInsert(t, h, gen->inter.uel);
            }
          });

        if(!force && UnsrtELSize(gen->inter.uel) == EDGECOUNT(nwp)){ // There are no ties in the initial network that are fixed.
          UnsrtELDestroy(gen->inter.uel);
          gen->inter.uel = NULL; // "Signal" that there is no discordance network.
        }else{
          gen->intertype = UnsrtELDyadGen;
          AddOnWtNetworkEdgeChange(nwp, (OnWtNetworkEdgeChange) WtDyadGenUpdate, gen, INT_MAX);
        }
      }
      break;
    default:
      error("Undefined dyad generator type.");
    }
}

void DyadGenUpgradeIntersect(DyadGen *gen){
  if(gen->intertype == UnsrtELDyadGen){
    gen->inter.hel = UnsrtELIntoHashEL(gen->inter.uel);
    gen->intertype = HashELDyadGen;
  }
}

DyadGen *DyadGenInitialize(DyadGenType type, void *dyads, void *track_nwp){
  DyadGen *gen = R_Calloc(1, DyadGen);
  gen->type = type;
  gen->intertype = NoELDyadGen;
  gen->sleeping = FALSE;
  gen->careless = TRUE;
  gen->inter.uel = NULL;

  switch(gen->type){
  case RandDyadGen:
    gen->nwp.b = dyads;
    gen->ndyads = DYADCOUNT(gen->nwp.b);
    break;
  case WtRandDyadGen:
    gen->nwp.w = dyads;
    gen->ndyads = DYADCOUNT(gen->nwp.w);
    break;
  case RLEBDM1DGen:
  case WtRLEBDM1DGen:
    // dyads must be a double ** that is shifted.
    gen->dyads.rlebdm = unpack_RLEBDM1D(dyads);
    gen->ndyads = gen->dyads.rlebdm.ndyads;
    break;
  case EdgeListGen:
  case WtEdgeListGen:
    // dyads must be an int ** that is shifted.
    gen->dyads.el = *(int **)dyads;
    gen->ndyads = *gen->dyads.el;
    (*(int **)dyads) += gen->ndyads*2 + 1;
    break;
  default:
    error("Undefined dyad generator type.");
  }

  if(track_nwp) DyadGenSetUpIntersect(gen, track_nwp, FALSE);

  for(unsigned int i = 0; i < kv_size(dyadgen_init_callbacks); i++)
    kv_A(dyadgen_init_callbacks, i).callback(gen, kv_A(dyadgen_init_callbacks, i).payload);

  return gen;
}


DyadGen *DyadGenInitializeR(SEXP pR, void *any_nwp, Rboolean el){
  // If there is a dyadgen element, use that; otherwise assume we're already getting a dyadgen list.
  SEXP dgR = getListElement(pR, "dyadgen");
  if(isNULL(dgR)) dgR = pR;

  DyadGenType type = asInteger(getListElement(dgR, "type"));

  void *track = el ? any_nwp : NULL;

  switch(type){
  case RandDyadGen:
  case WtRandDyadGen:
    return DyadGenInitialize(type, any_nwp, track);
  case RLEBDM1DGen:
  case WtRLEBDM1DGen:
    // RLEBDM1D's unpacking function expects a double **.
    {
      double *tmp = REAL(getListElement(dgR, "dyads"));
      return DyadGenInitialize(type, &tmp, track);
    }
  case EdgeListGen:
  case WtEdgeListGen:
    // RLEBDM1D's unpacking function expects an int **.
    {
      int *tmp = INTEGER(getListElement(dgR, "dyads"));
      return DyadGenInitialize(type, &tmp, track);
    }
  default:
    error("Undefined dyad generator type.");
  }
}


void DyadGenDestroy(DyadGen *gen){
  if(gen->intertype != NoELDyadGen){
    switch(gen->intertype){
    case UnsrtELDyadGen: UnsrtELDestroy(gen->inter.uel); break;
    case HashELDyadGen: HashELDestroy(gen->inter.hel); break;
    case NoELDyadGen: break;
    }

    switch(gen->type){
    case RLEBDM1DGen:
    case EdgeListGen:
      DeleteOnNetworkEdgeChange(gen->nwp.b, (OnNetworkEdgeChange) DyadGenUpdate, gen);
      break;
    case WtRLEBDM1DGen:
    case WtEdgeListGen:
      DeleteOnWtNetworkEdgeChange(gen->nwp.w, (OnWtNetworkEdgeChange) WtDyadGenUpdate, gen);
      break;
    case RandDyadGen:
    case WtRandDyadGen:
      break;
    default:
      error("Undefined dyad generator type.");
    }
  }
  R_Free(gen);
}


void DyadGenUpdate(Vertex tail, Vertex head, void *payload, Network *nwp, Rboolean edgestate){
  DyadGen *gen = payload;
  if(gen->sleeping) return;

  switch(gen->intertype){
  case UnsrtELDyadGen:
    if(gen->careless || DyadGenSearch(tail, head, gen)){
      if(edgestate) UnsrtELDelete(tail, head, gen->inter.uel); // Deleting
      else UnsrtELInsert(tail, head, gen->inter.uel); // Inserting
      if(gen->inter.uel->linsearch >= DYADGEN_MISSES_BEFORE_UPGRADE) DyadGenUpgradeIntersect(gen);
    }
    break;

  case HashELDyadGen:
    if(DyadGenSearch(tail, head, gen)){
      if(edgestate) HashELDelete(tail, head, gen->inter.hel); // Deleting
      else HashELInsert(tail, head, gen->inter.hel); // Inserting
    }
    break;

  case NoELDyadGen: break;
  }
}


void WtDyadGenUpdate(Vertex tail, Vertex head, double weight, void *payload, WtNetwork *nwp, double edgestate){
  DyadGen *gen = payload;
  if(gen->sleeping) return;

  switch(gen->intertype){
  case UnsrtELDyadGen:
    if(gen->careless || DyadGenSearch(tail, head, gen)){
      if(edgestate!=0 && weight==0) UnsrtELDelete(tail, head, gen->inter.uel); // Deleting
      else if(edgestate==0 && weight!=0) UnsrtELInsert(tail, head, gen->inter.uel); // Inserting
      if(gen->inter.uel->linsearch >= DYADGEN_MISSES_BEFORE_UPGRADE) DyadGenUpgradeIntersect(gen);
    }
    break;

  case HashELDyadGen:
    if(gen->careless || DyadGenSearch(tail, head, gen)){
      if(edgestate!=0 && weight==0) HashELDelete(tail, head, gen->inter.hel); // Deleting
      else if(edgestate==0 && weight!=0) HashELInsert(tail, head, gen->inter.hel); // Inserting
    }
    break;

  case NoELDyadGen: break;
  }
}

void AddOnDyadGenInit(OnDyadGenInit callback, void *payload){
  kv_push(OnDyadGenInitInfo, dyadgen_init_callbacks, ((OnDyadGenInitInfo) {.callback = callback, .payload = payload}));
}
void DeleteOnDyadGenInit(OnDyadGenInit callback, void *payload){
  unsigned int i = 0;
  while(i < kv_size(dyadgen_init_callbacks) &&
        (kv_A(dyadgen_init_callbacks, i).callback != callback ||
         kv_A(dyadgen_init_callbacks, i).payload != payload)) i++;

  if(i == kv_size(dyadgen_init_callbacks)) error("Attempting to delete a nonexistent DyadGen initialization callback.");

  kv_del_plug(dyadgen_init_callbacks, i);
}
