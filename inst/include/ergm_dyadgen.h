/*  File inst/include/ergm_dyadgen.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2024 Statnet Commons
 */
#ifndef _ERGM_DYADGEN_H_
#define _ERGM_DYADGEN_H_

#define STRICT_Wt_HEADERS
#include "ergm_edgetree.h"
#include "ergm_wtedgetree.h"
#include "ergm_rlebdm.h"
#include "ergm_edgelist.h"
#include "ergm_hash_edgelist.h"
#include "ergm_changestat.h"
#include "ergm_wtchangestat.h"

typedef enum{RandDyadGen, WtRandDyadGen, RLEBDM1DGen, WtRLEBDM1DGen, EdgeListGen, WtEdgeListGen} DyadGenType;
typedef enum{NoELDyadGen, UnsrtELDyadGen, HashELDyadGen} DyadGenInterType;
#define DYADGEN_MISSES_BEFORE_UPGRADE 8

typedef struct {
  DyadGenType type;
  DyadGenInterType intertype;
  union {
    Network *b;
    WtNetwork *w;
  } nwp;
  union {
    RLEBDM1D rlebdm;
    int *el;
  } dyads;
  Dyad ndyads;
  union {
    UnsrtEL *uel;
    HashEL *hel;
  } inter;
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

  if(gen->intertype == UnsrtELDyadGen){ /* If we are maintaining an unsorted edgelist (which also implies that we are *not* in a *RandDyadGen mode)... */
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
    if(extant) UnsrtELGetRand(tail, head, gen->inter.uel);
  }
}


static inline Edge DyadGenEdgecount(DyadGen *gen){
  switch(gen->intertype){
  case NoELDyadGen:
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
  case UnsrtELDyadGen: return UnsrtELSize(gen->inter.uel);
  case HashELDyadGen: return HashELSize(gen->inter.hel);
  }
  return 0; // Fix a warning.
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
      switch(gen->intertype){
      case NoELDyadGen: GetRandEdge(tail, head, gen->nwp.b); break;
      case UnsrtELDyadGen: UnsrtELGetRand(tail, head, gen->inter.uel); break;
      case HashELDyadGen: HashELGetRand(tail, head, gen->inter.hel); break;
      }
    }
    break;
  case WtRLEBDM1DGen:
  case WtEdgeListGen:
    {
      double dummy;
      switch(gen->intertype){
      case NoELDyadGen: WtGetRandEdge(tail, head, &dummy, gen->nwp.w); break;
      case UnsrtELDyadGen: UnsrtELGetRand(tail, head, gen->inter.uel); break;
      case HashELDyadGen: HashELGetRand(tail, head, gen->inter.hel); break;
      }
    }
    break;
  default:
    error("Undefined dyad generator type.");
  }
}


static inline void DyadGenRandNonedge(Vertex *tail, Vertex *head, DyadGen *gen){
  switch(gen->type){
  case RandDyadGen:
    GetRandNonedge(tail, head, gen->nwp.b);
    break;
  case WtRandDyadGen:
    WtGetRandNonedge(tail, head, gen->nwp.w);
    break;
  default:
    {
      Rboolean valued = gen->type == WtRLEBDM1DGen || gen->type == WtEdgeListGen;
      do{
        switch(gen->type){
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
      }while(valued ? WtEdgetreeSearch(*tail, *head, gen->nwp.w->outedges) : EdgetreeSearch(*tail, *head, gen->nwp.b->outedges));
    }
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
      switch(gen->intertype){
      case NoELDyadGen: GetRandEdge(tail, head, gen->nwp.b); break;
      case UnsrtELDyadGen: UnsrtELGetRand(tail, head, gen->inter.uel); break;
      case HashELDyadGen: HashELGetRand(tail, head, gen->inter.hel); break;
      }
    *weight = 1;
    break;
  case WtRLEBDM1DGen:
  case WtEdgeListGen:
      switch(gen->intertype){
      case NoELDyadGen: WtGetRandEdge(tail, head, weight, gen->nwp.w); break;
      case UnsrtELDyadGen:
        UnsrtELGetRand(tail, head, gen->inter.uel);
        *weight = WtGetEdge(*tail, *head, gen->nwp.w);
        break;
      case HashELDyadGen:
        HashELGetRand(tail, head, gen->inter.hel);
        *weight = WtGetEdge(*tail, *head, gen->nwp.w);
        break;
      }
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
  if(gen->intertype == NoELDyadGen) DyadGenSetUpIntersect(gen, NULL, TRUE);
}

static inline void DyadGenWake(DyadGen *gen){
  gen->sleeping = FALSE;
}

#endif
