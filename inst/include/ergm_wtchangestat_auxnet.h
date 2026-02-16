/*  File inst/include/ergm_wtchangestat_auxnet.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */
#ifndef _ERGM_WTCHANGESTAT_AUXNET_H_
#define _ERGM_WTCHANGESTAT_AUXNET_H_

#include "ergm_wtedgetree.h"
#include "ergm_wtchangestat.h"
#include "ergm_storage.h"
#include "ergm_wtmodel.h"
#include "ergm_wtchangestat_operator.h"
#include "ergm_edgelist.h"

/* Brief API description:

   The struct StoreAuxnet comprises the following information:

   * A pointer to the output network ownp.

   * A pointer to the input network inwp.

   * A pointer to the ModelTerm of the auxiliary that manages onwp.

   In an auxiliary, in addition to its usual i_, u_, and f_ functions,
   should provide:

   * a static inline function constructed by MAP_WtTOGGLE_FN() whose
     name is map_toggle_NAME() that computes the consequences of a
     toggle: return the number of induced toggles for the (tail,head)
     toggles (and 0 if none occurred), and put the tails and heads of
     those toggles in the respective arguments.

   * a macro (not a function or a variable!) named map_toggle_maxtoggles_NAME set to the maximum
     number of toggles produced by one input toggle.
*/

#define MAP_WtTOGGLE_FN(a) static inline unsigned int (a) (Vertex tail, Vertex head, double weight, double edgestate, struct StoreWtAuxnet_s *auxnet, Vertex *tails, Vertex *heads, double *weights)

typedef struct StoreWtAuxnet_s{WtNetwork *inwp, *onwp;
  WtModelTerm *mtp;
} StoreWtAuxnet;

#define MAP_WtTOGGLE(name, tail, head, weight, edgestate, auxnet, tails, heads, weights) map_toggle_ ## name(tail, head, weight, edgestate, auxnet, tails, heads, weights)

#define MAP_WtTOGGLE_THEN(name, tail, head, weight, edgestate, auxnet, tails, heads, weights) if(MAP_WtTOGGLE(name, tail, head, weight, edgestate, auxnet, tails, heads, weights))

#define I_WtAUXNET(init_onwp)                                   \
  ALLOC_AUX_STORAGE(1, StoreWtAuxnet, auxnet);			\
  auxnet->onwp = init_onwp;					\
  auxnet->inwp = nwp;						\
  auxnet->mtp = mtp;

#define MAP_WtTOGGLE_PROPAGATE *tails = tail; *heads = head; *weights = weight; return 1;
#define MAP_WtTOGGLE_PROPAGATE_IF(cond) if(cond){*tails = tail; *heads = head; *weights = weight; return 1;}else{return 0;}

#define ON_WtAUXNET(name)                                               \
  WtI_CHANGESTAT_FN(i_on ## name){                                      \
    GET_AUX_STORAGE(StoreWtAuxnet, auxnet);                             \
    WtModel *m = STORAGE = WtModelInitialize(getListElement(mtp->R, "submodel"),  NULL, auxnet->onwp, FALSE); \
    WtDELETE_IF_UNUSED_IN_SUBMODEL(u_func, m);                          \
    WtDELETE_IF_UNUSED_IN_SUBMODEL(z_func, m);                          \
    /* SELECT_C_OR_D_BASED_ON_SUBMODEL(m); */                           \
  }                                                                     \
                                                                        \
  WtC_CHANGESTAT_FN(c_on ## name){                                      \
    GET_STORAGE(WtModel, m);                                            \
    GET_AUX_STORAGE(StoreWtAuxnet, auxnet);                             \
                                                                        \
    Vertex tails[map_toggle_maxtoggles_ ## name], heads[map_toggle_maxtoggles_ ## name]; \
    double weights[map_toggle_maxtoggles_ ## name];                     \
    if(map_toggle_maxtoggles_ ## name == 1){ /* One of these should get optimized away by the compiler. */ \
      MAP_WtTOGGLE_THEN(name, tail, head, weight, edgestate, auxnet, tails, heads, weights){ \
        double *tmp = m->workspace;                                     \
        m->workspace = CHANGE_STAT;                                     \
        WtChangeStats1(*tails, *heads, *weights, auxnet->onwp, m, WtGETWT(*tails, *heads, auxnet->onwp)); \
        m->workspace = tmp;                                             \
      }                                                                 \
    }else{                                                              \
      unsigned int ntoggles = MAP_WtTOGGLE(name, tail, head, weight, edgestate, auxnet, tails, heads, weights); \
      if(ntoggles){                                                     \
        double *tmp = m->workspace;                                     \
        m->workspace = CHANGE_STAT;                                     \
        WtChangeStats(ntoggles, tails, heads, weights, auxnet->onwp, m); \
        m->workspace = tmp;                                             \
      }                                                                 \
    }                                                                   \
  }                                                                     \
                                                                        \
  WtZ_CHANGESTAT_FN(z_on ## name){                                      \
    GET_STORAGE(WtModel, m);                                            \
                                                                        \
    double *tmp = m->workspace;                                         \
    m->workspace = CHANGE_STAT;                                         \
    WtZStats(nwp, m, skip_s);                                           \
    m->workspace = tmp;                                                 \
  }                                                                     \
                                                                        \
  WtF_CHANGESTAT_FN(f_on ## name){                                      \
    GET_STORAGE(WtModel, m);                                            \
    GET_AUX_STORAGE(StoreWtAuxnet, auxnet);                             \
    WtModelDestroy(auxnet->onwp, m);                                    \
    STORAGE = NULL;                                                     \
  }

#endif // _ERGM_CHANGESTAT_AUXNET_H_
