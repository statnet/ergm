/*  File inst/include/ergm_dyadgen.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2022 Statnet Commons
 */
#ifndef _ERGM_DYADGEN_H_
#define _ERGM_DYADGEN_H_

#define STRICT_Wt_HEADERS
#include "ergm_edgetree.h"
#include "ergm_wtedgetree.h"
#include "ergm_rlebdm.h"
#include "ergm_edgelist.h"
#include "ergm_unsorted_edgelist.h"
#include "ergm_changestat.h"
#include "ergm_wtchangestat.h"

typedef enum{RandDyadGen, WtRandDyadGen, RLEBDM1DGen, WtRLEBDM1DGen, EdgeListGen, WtEdgeListGen} DyadGenType;

typedef struct {
  DyadGenType type;
  union {
    Network *b;
    WtNetwork *w;
  } nwp;
  union {
    RLEBDM1D rlebdm;
    int *el;
  } dyads;
  Dyad ndyads;
  UnsrtEL *intersect;
  Rboolean sleeping;
} DyadGen;

void DyadGenSetUpIntersect(DyadGen *gen, void *track_nwp, Rboolean force);
DyadGen *DyadGenInitialize(DyadGenType type, void *dyads, void *track_nwp);
DyadGen *DyadGenInitializeR(SEXP pR, void *any_nwp, Rboolean el);
void DyadGenDestroy(DyadGen *gen);

void DyadGenUpdate(Vertex tail, Vertex head, DyadGen *gen, Network *nwp, Rboolean edgestate);
void WtDyadGenUpdate(Vertex tail, Vertex head, double weight, DyadGen *gen, WtNetwork *nwp, double edgestate);

static inline void DyadGenRandDyad(Vertex *tail, Vertex *head, DyadGen *gen){
  switch(gen->type){
  case RandDyadGen:
    GetRandDyad(tail, head, gen->nwp.b);
    break;
  case WtRandDyadGen:
    GetRandDyad(tail, head, gen->nwp.w);
    break;
  case RLEBDM1DGen:
  case WtRLEBDM1DGen:
    GetRandRLEBDM1D(tail, head, &gen->dyads.rlebdm);
    break;
  case EdgeListGen:
  case WtEdgeListGen:
    {
      int ndyads = gen->ndyads;
      int *list = gen->dyads.el + 1;
      Edge rane = unif_rand() * ndyads;
      *tail = list[rane];
      *head = list[ndyads+rane];
    }
    break;
  default:
    error("Undefined dyad generator type.");
  }

  if(gen->intersect){ /* If we are maintaining an unsorted edgelist (which also implies that we are *not* in a *RandDyadGen mode)... */
    /* Use the appropriate function to check if we had selected an extant edge, */
    Rboolean extant;
    switch(gen->type){
    case RLEBDM1DGen:
    case EdgeListGen:
      extant = EdgetreeSearch(*tail, *head, gen->nwp.b->outedges)!=0;
      break;
    case WtRLEBDM1DGen:
    case WtEdgeListGen:
      extant = WtEdgetreeSearch(*tail, *head, gen->nwp.w->outedges)!=0;
      break;
    default:
      error("Undefined dyad generator type.");
    }

    /* ... and if so, reselect it from the unsorted edgelist. */
    if(extant) UnsrtELGetRand(tail, head, gen->intersect);
  }
}


static inline Edge DyadGenEdgecount(DyadGen *gen){
  if(gen->intersect){
    return gen->intersect->nedges;
  }else{
    switch(gen->type){
    case RandDyadGen:
    case RLEBDM1DGen:
    case EdgeListGen:
      return EDGECOUNT(gen->nwp.b);
    case WtRandDyadGen:
    case WtRLEBDM1DGen:
    case WtEdgeListGen:
      return EDGECOUNT(gen->nwp.w);
    default:
      error("Undefined dyad generator type.");
    }
  }
}


static inline void DyadGenRandEdge(Vertex *tail, Vertex *head, DyadGen *gen){
  switch(gen->type){
  case RandDyadGen:
    GetRandEdge(tail, head, gen->nwp.b);
    break;
  case WtRandDyadGen:
    {
      double dummy;
      WtGetRandEdge(tail, head, &dummy, gen->nwp.w);
    }
    break;
  case RLEBDM1DGen:
  case EdgeListGen:
    {
      if(gen->intersect) UnsrtELGetRand(tail, head, gen->intersect);
      else GetRandEdge(tail, head, gen->nwp.b);
    }
    break;
  case WtRLEBDM1DGen:
  case WtEdgeListGen:
    {
      double dummy;
      if(gen->intersect) UnsrtELGetRand(tail, head, gen->intersect);
      else WtGetRandEdge(tail, head, &dummy, gen->nwp.w);
    }
    break;
  default:
    error("Undefined dyad generator type.");
  }
}


static inline void DyadGenRandWtEdge(Vertex *tail, Vertex *head, double *weight, DyadGen *gen){
  switch(gen->type){
  case RandDyadGen:
    GetRandEdge(tail, head, gen->nwp.b);
    *weight = 1;
    break;
  case WtRandDyadGen:
    WtGetRandEdge(tail, head, weight, gen->nwp.w);
    break;
  case RLEBDM1DGen:
  case EdgeListGen:
    if(gen->intersect) UnsrtELGetRand(tail, head, gen->intersect);
    else GetRandEdge(tail, head, gen->nwp.b);
    *weight = 1;
    break;
  case WtRLEBDM1DGen:
  case WtEdgeListGen:
    if(gen->intersect){
      UnsrtELGetRand(tail, head, gen->intersect);
      *weight = WtGetEdge(*tail, *head, gen->nwp.w);
    }
    else WtGetRandEdge(tail, head, weight, gen->nwp.w);
    break;
  default:
    error("Undefined dyad generator type.");
  }
}


static inline Rboolean DyadGenSearch(Vertex tail, Vertex head, DyadGen *gen){
  switch(gen->type){
  case RandDyadGen:
  case WtRandDyadGen:
    return TRUE;
  case RLEBDM1DGen:
  case WtRLEBDM1DGen:
    return GetRLEBDM1D(tail, head, &gen->dyads.rlebdm);
  case EdgeListGen:
  case WtEdgeListGen:
    return iEdgeListSearch(tail, head, gen->dyads.el);
  default:
    error("Undefined dyad generator type.");
  }
}

static inline void DyadGenSleep(DyadGen *gen){
  gen->sleeping = TRUE;

  /* If this DyadGen is asked to sleep, that means that someone else
     might add edges that don't "belong" to this DyadGen, which means
     that it must keep track of what they are.

     FIXME: Maybe we should *always* maintain an UnsrtEL? Need to benchmark. */
  if(!gen->intersect) DyadGenSetUpIntersect(gen, NULL, TRUE);
}

static inline void DyadGenWake(DyadGen *gen){
  gen->sleeping = FALSE;
}

#endif
