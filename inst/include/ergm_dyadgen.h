#ifndef _ERGM_DYADGEN_H_
#define _ERGM_DYADGEN_H_

#define STRICT_WT_HEADERS
#include "ergm_edgetree.h"
#include "ergm_wtedgetree.h"
#include "ergm_rlebdm.h"
#include "ergm_edgelist.h"

enum DyadGenType {RandDyadGen, WtRandDyadGen, RLEBDM1DGen, EdgeListGen};

typedef struct {
  enum DyadGenType type;
  union {
    Network *nwp;
    WtNetwork *wtnwp;
    RLEBDM1D rlebdm;
    int *el;
  } data;
  Dyad ndyads;
} DyadGen;

static inline void GenRandDyad(Vertex *tail, Vertex *head, DyadGen *gen){
  switch(gen->type){
  case RandDyadGen:
    GetRandDyad(tail, head, gen->data.nwp);
    break;
  case WtRandDyadGen:
    GetRandDyad(tail, head, gen->data.wtnwp);
    break;
  case RLEBDM1DGen:
    GetRandRLEBDM1D(tail, head, &gen->data.rlebdm);
    break;
  case EdgeListGen:
    {
      int ndyads = gen->ndyads;
      int *list = gen->data.el + 1;
      Edge rane = unif_rand() * ndyads;
      *tail = list[rane];
      *head = list[ndyads+rane];
    }
    break;
  default:
    error("Undefined dyad generator type.");
  }
}

static inline Rboolean GetDyadGen(Vertex tail, Vertex head, DyadGen *gen){
  switch(gen->type){
  case RandDyadGen:
  case WtRandDyadGen:
    return TRUE;
    break;
  case RLEBDM1DGen:
    return GetRLEBDM1D(tail, head, &gen->data.rlebdm);
    break;
  case EdgeListGen:
    return iEdgeListSearch(tail, head, gen->data.el);
    break;
  default:
    error("Undefined dyad generator type.");
  }
}

DyadGen *DyadGenInitialize(enum DyadGenType type, void *data);
DyadGen *DyadGenInitializeR(SEXP pR, void *any_nwp);
void DyadGenDestroy(DyadGen *gen);

#endif
