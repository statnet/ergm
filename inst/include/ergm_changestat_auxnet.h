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

/* Brief API description:

   The struct StoreAuxnet comprises the following information:
   
   * A pointer to the output network ownp.

   * A pointer to the input network inwp.

   * A pointer to the ModelTerm of the auxiliary that manages onwp.

   * A pointer to a function that calculates the consequences of a
     particular toggle.

   The map_toggle() function takes toggle tail, head, a StoreAuxnet, a
   pointer to an unsigned int, and a pointer to a vector of tails and
   heads, and should heave as follows:

   * If tail == 0, set *ntoggles to the maximum number of toggles in
     the output network that a single toggle in the input network
     could induce.

   * Otherwise, set *ntoggles to the number of induced toggles for the
     (tail,head) toggles (and 0 if none occurred), and put the tails
     and heads of those toggles in the respective arguments.

   It is generally recommended to call it via the MAP_TOGGLE macro or
   MAP_TOGGLE_1 if the maximum number of induced toggles is known a
   priori to be 1.
*/

#define MAP_TOGGLE_ARGS (Vertex tail, Vertex head, struct StoreAuxnet_s *auxnet, unsigned int *ntoggles, Vertex *tails, Vertex *heads)
#define MAP_TOGGLE_FN(a) static void (a) MAP_TOGGLE_ARGS

typedef struct StoreAuxnet_s{Network *inwp, *onwp;
  ModelTerm *mtp;
  void (*map_toggle) MAP_TOGGLE_ARGS;
} StoreAuxnet;

#define MAP_TOGGLE(tail, head, auxnet, ntoggles, tails, heads) auxnet->map_toggle(tail, head, auxnet, &ntoggles, tails, heads)

#define MAP_TOGGLE_1(tail, head, auxnet, ntoggles, tails, heads)	\
  unsigned int ntoggles;						\
  Vertex tails[1], heads[1];						\
  MAP_TOGGLE(tail, head, auxnet, ntoggles, tails, heads);

#define I_AUXNET(init_onwp, init_map_toggle)			\
  ALLOC_AUX_STORAGE(1, StoreAuxnet, auxnet);			\
  auxnet->onwp = init_onwp;					\
  auxnet->inwp = nwp;						\
  auxnet->map_toggle = init_map_toggle;				\
  auxnet->mtp = mtp;

#define MAP_TOGGLE_PROPAGATE *ntoggles = 1; *tails = tail; *heads = head; return;
#define MAP_TOGGLE_PROPAGATE_IF(cond) if(cond){*ntoggles = 1; *tails = tail; *heads = head;}else{*ntoggles = 0;}; return;

#define MAP_TOGGLE_MAXTOGGLES(maxtoggles) if(tail==0){*ntoggles = maxtoggles;return;}



#endif // _ERGM_CHANGESTAT_AUXNET_H_
