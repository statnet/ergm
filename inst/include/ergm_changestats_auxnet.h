/*  File inst/include/ergm_changestats_auxnet.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2021 Statnet Commons
 */
#ifndef _ERGM_CHANGESTATS_AUXNET_H_
#define _ERGM_CHANGESTATS_AUXNET_H_

#include "ergm_changestat_auxnet.h"
#include "ergm_dyad_hashmap.h"

typedef struct{StoreDyadSet *nwp; int *ref_el;} StoreDyadSetAndRefEL;

#define map_toggle_maxtoggles__discord_net_Network 1
MAP_TOGGLE_FN(map_toggle__discord_net_Network){
  MAP_TOGGLE_PROPAGATE;
}

#define map_toggle_maxtoggles__intersect_net_Network 1
MAP_TOGGLE_FN(map_toggle__intersect_net_Network){
  ModelTerm *mtp = auxnet->mtp;
  int *ref_el = IINPUT_PARAM;
  MAP_TOGGLE_PROPAGATE_IF(iEdgeListSearch(tail, head, ref_el));
}

#define map_toggle_maxtoggles__union_net_Network 1
MAP_TOGGLE_FN(map_toggle__union_net_Network){
  ModelTerm *mtp = auxnet->mtp;
  int *ref_el = IINPUT_PARAM;
  MAP_TOGGLE_PROPAGATE_IF(!iEdgeListSearch(tail, head, ref_el));
}

#define map_toggle_maxtoggles__blockdiag_net 1
MAP_TOGGLE_FN(map_toggle__blockdiag_net){
  ModelTerm *mtp = auxnet->mtp;
  int *b = IINPUT_PARAM-1; // tail and head are indexed from 1.
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

#define map_toggle_maxtoggles__undir_net 1
MAP_TOGGLE_FN(map_toggle__undir_net){
  ModelTerm *mtp = auxnet->mtp;
  unsigned int rule = IINPUT_PARAM[0];
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

#define map_toggle_maxtoggles__filter_formula_net 1
MAP_TOGGLE_FN(map_toggle__filter_formula_net){
  ModelTerm *mtp = auxnet->mtp;
  GET_STORAGE(Model, m);
  ChangeStats1(tail, head, auxnet->inwp, m, edgestate);
  MAP_TOGGLE_PROPAGATE_IF(*(m->workspace)!=0);
}

#define map_toggle_maxtoggles__subgraph_net 1
MAP_TOGGLE_FN(map_toggle__subgraph_net){
  Network *nwp = auxnet->inwp;
  ModelTerm *mtp = auxnet->mtp;
  GET_STORAGE(int*, thmap);

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
