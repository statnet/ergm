/*  File src/changestats_test.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2017 Statnet Commons
 */
#ifndef _ERGM_CHANGESTAT_AUXNET_H_
#define _ERGM_CHANGESTAT_AUXNET_H_

#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_storage.h"
#include "ergm_changestat_operator.h"

/* Brief API description:

   The struct StoreAuxnet comprises the following information:
   
   * A pointer to the output network ownp.

   * A pointer to the input network inwp.

   * A pointer to the ModelTerm of the auxiliary that manages onwp.

   In addition, a pair of static inline functions constructed by
   MAP_TOGGLE_FN() and by MAP_TOGGLE_MAXTOGGLES_FN() should be
   provided in a header file, with their names being map_toggle_NAME
   and map_toggle_maxtoggles_NAME that compute consequences of a
   toggle and the maximum number of toggles produced by one input
   toggle as follows:

   * If tail == 0, set return the maximum number of toggles in the
     output network that a single toggle in the input network could
     induce.

   * Otherwise, return the number of induced toggles for the
     (tail,head) toggles (and 0 if none occurred), and put the tails
     and heads of those toggles in the respective arguments.

   It is generally recommended to call it via the MAP_TOGGLE macro or
   MAP_TOGGLE_1_THEN if the maximum number of induced toggles is known
   a priori to be 1.
*/

#define MAP_TOGGLE_FN(a) static inline unsigned int (a) (Vertex tail, Vertex head, Rboolean edgeflag, struct StoreAuxnet_s *auxnet, Vertex *tails, Vertex *heads)
#define MAP_TOGGLE_MAXTOGGLES_FN(a) static inline unsigned int (a) (struct StoreAuxnet_s *auxnet)

typedef struct StoreAuxnet_s{Network *inwp, *onwp;
  ModelTerm *mtp;
} StoreAuxnet;

#define MAP_TOGGLE(name, tail, head, edgeflag, auxnet, tails, heads) map_toggle_ ## name(tail, head, edgeflag, auxnet, tails, heads)

#define MAP_TOGGLE_THEN(name, tail, head, edgeflag, auxnet, tails, heads) if(MAP_TOGGLE(name, tail, head, edgeflag, auxnet, tails, heads))

#define MAP_TOGGLE_1_THEN(name, tail, head, edgeflag, auxnet, tails, heads) \
  Vertex tails[1], heads[1];                                          \
  MAP_TOGGLE_THEN(name, tail, head, edgeflag, auxnet, tails, heads)

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
    double *inputs = INPUT_PARAM + 1;                                   \
    STORAGE = unpack_Model_as_double(&inputs, auxnet->onwp);            \
  }                                                                     \
                                                                        \
  C_CHANGESTAT_FN(c_on ## name){                                        \
    GET_STORAGE(Model, m);                                              \
    GET_AUX_STORAGE(StoreAuxnet, auxnet);                               \
                                                                        \
    MAP_TOGGLE_1_THEN(name, tail, head, edgeflag, auxnet, tails, heads){ \
      ChangeStats(1, tails, heads, auxnet->onwp, m);                    \
      memcpy(CHANGE_STAT, m->workspace, N_CHANGE_STATS*sizeof(double)); \
    }                                                                   \
  }                                                                     \
                                                                        \
  U_CHANGESTAT_FN(u_on ## name){                                        \
    GET_STORAGE(Model, m);                                              \
    GET_AUX_STORAGE(StoreAuxnet, auxnet);                               \
                                                                        \
    MAP_TOGGLE_1_THEN(name, tail, head, edgeflag, auxnet, tails, heads) UPDATE_STORAGE(*tails, *heads, auxnet->onwp, m, NULL, edgeflag); \
  }                                                                     \
                                                                        \
  F_CHANGESTAT_FN(f_on ## name){                                        \
    GET_STORAGE(Model, m);                                              \
    GET_AUX_STORAGE(StoreAuxnet, auxnet);                               \
    ModelDestroy(auxnet->onwp, m);                                      \
    STORAGE = NULL;                                                     \
  }

#endif // _ERGM_CHANGESTAT_AUXNET_H_
