/*  File inst/include/ergm_changestat_auxnet.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2024 Statnet Commons
 */
#ifndef _ERGM_CHANGESTAT_AUXNET_H_
#define _ERGM_CHANGESTAT_AUXNET_H_

#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_storage.h"
#include "ergm_model.h"
#include "ergm_changestat_operator.h"
#include "ergm_edgelist.h"

/* Brief API description:

   The struct StoreAuxnet comprises the following information:
   
   * A pointer to the output network ownp.

   * A pointer to the input network inwp.

   * A pointer to the ModelTerm of the auxiliary that manages onwp.

   In an auxiliary, in addition to its usual i_, u_, and f_ functions,
   should provide:

   * a static inline function constructed by MAP_TOGGLE_FN() whose
     name is map_toggle_NAME() that computes the consequences of a
     toggle: return the number of induced toggles for the (tail,head)
     toggles (and 0 if none occurred), and put the tails and heads of
     those toggles in the respective arguments.

   * a macro (not a function or a variable!) named map_toggle_maxtoggles_NAME set to the maximum
     number of toggles produced by one input toggle.
*/

#define MAP_TOGGLE_FN(a) static inline unsigned int (a) (Vertex tail, Vertex head, Rboolean edgestate, struct StoreAuxnet_s *auxnet, Vertex *tails, Vertex *heads)

typedef struct StoreAuxnet_s{Network *inwp, *onwp;
  ModelTerm *mtp;
} StoreAuxnet;

#define MAP_TOGGLE(name, tail, head, edgestate, auxnet, tails, heads) map_toggle_ ## name(tail, head, edgestate, auxnet, tails, heads)

#define MAP_TOGGLE_THEN(name, tail, head, edgestate, auxnet, tails, heads) if(MAP_TOGGLE(name, tail, head, edgestate, auxnet, tails, heads))

#define I_AUXNET(init_onwp)                                     \
  ALLOC_AUX_STORAGE(1, StoreAuxnet, auxnet);			\
  auxnet->onwp = init_onwp;					\
  auxnet->inwp = nwp;						\
  auxnet->mtp = mtp;

#define MAP_TOGGLE_PROPAGATE *tails = tail; *heads = head; return 1;
#define MAP_TOGGLE_PROPAGATE_IF(cond) if(cond){*tails = tail; *heads = head; return 1;}else{return 0;}

#define ON_AUXNET(name)                                                 \
  I_CHANGESTAT_FN(i_on ## name){                                        \
    GET_AUX_STORAGE(StoreAuxnet, auxnet);                               \
    Model *m = STORAGE = ModelInitialize(getListElement(mtp->R, "submodel"),  NULL, auxnet->onwp, FALSE); \
    DELETE_IF_UNUSED_IN_SUBMODEL(u_func, m);                            \
    DELETE_IF_UNUSED_IN_SUBMODEL(z_func, m);                            \
    /* SELECT_C_OR_D_BASED_ON_SUBMODEL(m); */                           \
  }                                                                     \
                                                                        \
  C_CHANGESTAT_FN(c_on ## name){                                        \
    GET_STORAGE(Model, m);                                              \
    GET_AUX_STORAGE(StoreAuxnet, auxnet);                               \
                                                                        \
    Vertex tails[map_toggle_maxtoggles_ ## name], heads[map_toggle_maxtoggles_ ## name]; \
    if(map_toggle_maxtoggles_ ## name == 1){ /* One of these should get optimized away by the compiler. */ \
      MAP_TOGGLE_THEN(name, tail, head, edgestate, auxnet, tails, heads){ \
        double *tmp = m->workspace;                                     \
        m->workspace = CHANGE_STAT;                                     \
        ChangeStats1(*tails, *heads, auxnet->onwp, m, IS_OUTEDGE(*tails, *heads, auxnet->onwp)); \
        m->workspace = tmp;                                             \
      }                                                                 \
    }else{                                                              \
      unsigned int ntoggles = MAP_TOGGLE(name, tail, head, edgestate, auxnet, tails, heads); \
      if(ntoggles){                                                     \
        double *tmp = m->workspace;                                     \
        m->workspace = CHANGE_STAT;                                     \
        ChangeStats(ntoggles, tails, heads, auxnet->onwp, m);           \
        m->workspace = tmp;                                             \
      }                                                                 \
    }                                                                   \
  }                                                                     \
                                                                        \
  Z_CHANGESTAT_FN(z_on ## name){                                        \
    GET_STORAGE(Model, m);                                              \
                                                                        \
    double *tmp = m->workspace;                                         \
    m->workspace = CHANGE_STAT;                                         \
    ZStats(nwp, m, skip_s);                                             \
    m->workspace = tmp;                                                 \
  }                                                                     \
                                                                        \
  F_CHANGESTAT_FN(f_on ## name){                                        \
    GET_STORAGE(Model, m);                                              \
    GET_AUX_STORAGE(StoreAuxnet, auxnet);                               \
    ModelDestroy(auxnet->onwp, m);                                      \
    STORAGE = NULL;                                                     \
  }

#endif // _ERGM_CHANGESTAT_AUXNET_H_
