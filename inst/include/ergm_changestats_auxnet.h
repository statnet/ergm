/*  File src/changestats_test.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2017 Statnet Commons
 */
#ifndef _ERGM_CHANGESTATS_AUXNET_H_
#define _ERGM_CHANGESTATS_AUXNET_H_

#include "ergm_changestat_auxnet.h"
#include "ergm_dyad_hashmap.h"

typedef struct{StoreDyadSet *nwp; double *ref_el;} StoreDyadSetAndRefEL;

I_CHANGESTAT_FN(i__isociomatrix);
U_CHANGESTAT_FN(u__isociomatrix);
F_CHANGESTAT_FN(f__isociomatrix);

I_CHANGESTAT_FN(i__discord_net_DyadSet);
U_CHANGESTAT_FN(u__discord_net_DyadSet);
F_CHANGESTAT_FN(f__discord_net_DyadSet);

I_CHANGESTAT_FN(i__intersect_net_DyadSet);
U_CHANGESTAT_FN(u__intersect_net_DyadSet);
F_CHANGESTAT_FN(f__intersect_net_DyadSet);

I_CHANGESTAT_FN(i__intersect_net_toggles_in_list_DyadSet);
U_CHANGESTAT_FN(u__intersect_net_toggles_in_list_DyadSet);
F_CHANGESTAT_FN(f__intersect_net_toggles_in_list_DyadSet);

I_CHANGESTAT_FN(i__union_net_DyadSet);
U_CHANGESTAT_FN(u__union_net_DyadSet);
F_CHANGESTAT_FN(f__union_net_DyadSet);


MAP_TOGGLE_MAXTOGGLES_FN(map_toggle_maxtoggles__discord_net_Network){
  return 1;
}

MAP_TOGGLE_FN(map_toggle__discord_net_Network){
  MAP_TOGGLE_PROPAGATE;
}

MAP_TOGGLE_MAXTOGGLES_FN(map_toggle_maxtoggles__intersect_net_Network){
  return 1;
}

MAP_TOGGLE_FN(map_toggle__intersect_net_Network){
  ModelTerm *mtp = auxnet->mtp;
  double *ref_el = INPUT_PARAM + 1;
  MAP_TOGGLE_PROPAGATE_IF(dEdgeListSearch(tail, head, ref_el));
}

MAP_TOGGLE_MAXTOGGLES_FN(map_toggle_maxtoggles__union_net_Network){
  return 1;
}

MAP_TOGGLE_FN(map_toggle__union_net_Network){
  ModelTerm *mtp = auxnet->mtp;
  double *ref_el = INPUT_PARAM + 1;
  MAP_TOGGLE_PROPAGATE_IF(!dEdgeListSearch(tail, head, ref_el));
}

MAP_TOGGLE_MAXTOGGLES_FN(map_toggle_maxtoggles__blockdiag_net){
  return 1;
}

MAP_TOGGLE_FN(map_toggle__blockdiag_net){
  ModelTerm *mtp = auxnet->mtp;
  double *b = INPUT_PARAM;
  MAP_TOGGLE_PROPAGATE_IF(b[tail]==b[head]);
}

#define __undir_net_totoggle				\
  Rboolean totoggle;					\
  switch(rule){						\
  case 1: /* weak */					\
    totoggle = !IS_OUTEDGE(head,tail);			\
    break;						\
  case 2: /* strong */					\
    totoggle = IS_OUTEDGE(head,tail);			\
    break;						\
  case 3: /* upper */					\
    totoggle = tail<=head;				\
    break;						\
  case 4: /* lower */					\
    totoggle = tail>=head;				\
    break;						\
  default: /* never reached, but avoids a warning */	\
    totoggle = FALSE;					\
  }

MAP_TOGGLE_MAXTOGGLES_FN(map_toggle_maxtoggles__undir_net){
  return 1;
}

MAP_TOGGLE_FN(map_toggle__undir_net){
  ModelTerm *mtp = auxnet->mtp;
  unsigned int rule = INPUT_PARAM[1];
  Network *nwp = auxnet->inwp;
  __undir_net_totoggle;
  if(totoggle){
    *tails = MIN(tail,head);
    *heads = MAX(tail,head);
    return 1;
  }else{
    return 0;
  }
}

MAP_TOGGLE_MAXTOGGLES_FN(map_toggle_maxtoggles__filter_formula_net){
  return 1;
}

MAP_TOGGLE_FN(map_toggle__filter_formula_net){
  ModelTerm *mtp = auxnet->mtp;
  GET_STORAGE(Model, m);
  ChangeStats(1, &tail, &head, auxnet->inwp, m);
  MAP_TOGGLE_PROPAGATE_IF(*(m->workspace)!=0);
}

MAP_TOGGLE_MAXTOGGLES_FN(map_toggle_maxtoggles__subgraph_net){
  return 1;
}

MAP_TOGGLE_FN(map_toggle__subgraph_net){
  Network *nwp = auxnet->inwp;
  ModelTerm *mtp = auxnet->mtp;
  GET_STORAGE(double*, thmap);

  Vertex st = thmap[0][tail];
  Vertex sh = thmap[1][head];
  if(!DIRECTED && (st==0 || sh==0)){
    st = thmap[0][head];
    sh = thmap[1][tail];
  }
  if(st!=0 && sh!=0){
    *tails = st;
    *heads = sh;
    return 1;
  }else{
    return 0;
  }
}


#endif
