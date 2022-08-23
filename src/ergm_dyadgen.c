/*  File src/ergm_dyadgen.c in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2022 Statnet Commons
 */
#include "ergm_dyadgen.h"
#include "ergm_Rutil.h"
#include "ergm_changestat.h"
#include "ergm_wtchangestat.h"

void DyadGenSetUpIntersect(DyadGen *gen, void *track_nwp, Rboolean force){
    switch(gen->type){
    case RandDyadGen:
    case WtRandDyadGen:
      gen->intersect = NULL;
      break;
    case EdgeListGen:
    case RLEBDM1DGen:
      {
        Network *nwp = track_nwp ? track_nwp : gen->nwp.b;
        gen->nwp.b = nwp;
        gen->intersect = UnsrtELInitialize(0, NULL, NULL, FALSE);
        EXEC_THROUGH_NET_EDGES(t, h, e, {
            if(DyadGenSearch(t, h, gen)){
              UnsrtELInsert(t, h, gen->intersect);
            }
          });
        
        if(!force && gen->intersect->nedges==EDGECOUNT(nwp)){ // There are no ties in the initial network that are fixed.
          UnsrtELDestroy(gen->intersect);
          gen->intersect = NULL; // "Signal" that there is no discordance network.
        }else{
          AddOnNetworkEdgeChange(nwp, (OnNetworkEdgeChange) DyadGenUpdate, gen, INT_MAX);
        }
      }
      break;
    case WtEdgeListGen:
    case WtRLEBDM1DGen:
      {
        WtNetwork *nwp = track_nwp ? track_nwp : gen->nwp.w;
        gen->nwp.w = nwp;
        gen->intersect = UnsrtELInitialize(0, NULL, NULL, FALSE);
        WtEXEC_THROUGH_NET_EDGES(t, h, w, e, {
            (void) e;
            if(w!=0 && DyadGenSearch(t, h, gen)){
              UnsrtELInsert(t, h, gen->intersect);
            }
          });
        
      if(!force && gen->intersect->nedges==EDGECOUNT(nwp)){ // There are no ties in the initial network that are fixed.
        UnsrtELDestroy(gen->intersect);
        gen->intersect = NULL; // "Signal" that there is no discordance network.
      }else{
        AddOnWtNetworkEdgeChange(nwp, (OnWtNetworkEdgeChange) WtDyadGenUpdate, gen, INT_MAX);
      }
      }
      break;
    default:
      error("Undefined dyad generator type.");
    }
}

DyadGen *DyadGenInitialize(DyadGenType type, void *dyads, void *track_nwp){
  DyadGen *gen = Calloc(1, DyadGen);
  gen->type = type;
  gen->sleeping = FALSE;

  switch(gen->type){
  case RandDyadGen:
    gen->nwp.b = dyads;
    gen->ndyads = DYADCOUNT(gen->nwp.b);
    gen->intersect = NULL;
    break;
  case WtRandDyadGen:
    gen->nwp.w = dyads;
    gen->ndyads = DYADCOUNT(gen->nwp.w);
    gen->intersect = NULL;
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
  if(gen->intersect){
    switch(gen->type){
    case RLEBDM1DGen:
    case EdgeListGen:
      if(gen->intersect){
        UnsrtELDestroy(gen->intersect);
        DeleteOnNetworkEdgeChange(gen->nwp.b, (OnNetworkEdgeChange) DyadGenUpdate, gen);
      }
      break;
    case WtRLEBDM1DGen:
    case WtEdgeListGen:
      if(gen->intersect){
        UnsrtELDestroy(gen->intersect);
        DeleteOnWtNetworkEdgeChange(gen->nwp.w, (OnWtNetworkEdgeChange) WtDyadGenUpdate, gen);
      }
      break;
    case RandDyadGen:
    case WtRandDyadGen:
      break;
    default:
      error("Undefined dyad generator type.");
    }
  }
  Free(gen);
}


void DyadGenUpdate(Vertex tail, Vertex head, DyadGen *gen, Network *nwp, Rboolean edgestate){
  if(gen->sleeping) return;
  if(edgestate) UnsrtELDelete(tail, head, gen->intersect); // Deleting
  else UnsrtELInsert(tail, head, gen->intersect); // Inserting
}


void WtDyadGenUpdate(Vertex tail, Vertex head, double weight, DyadGen *gen, WtNetwork *nwp, double edgestate){
  if(gen->sleeping) return;
  if(edgestate!=0 && weight==0) UnsrtELDelete(tail, head, gen->intersect); // Deleting
  else if(edgestate==0 && weight!=0) UnsrtELInsert(tail, head, gen->intersect); // Inserting
}
