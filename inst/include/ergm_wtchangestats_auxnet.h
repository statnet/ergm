/*  File inst/include/ergm_wtchangestats_auxnet.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#ifndef _ERGM_CHANGESTATS_AUXNET_H_
#define _ERGM_CHANGESTATS_AUXNET_H_

#include "ergm_wtchangestat_auxnet.h"
#include "ergm_dyad_hashmap.h"

typedef struct{StoreStrictDyadSet *nwp; int *ref_el;} StoreStrictDyadSetAndRefEL;

/* #define map_toggle_maxtoggles__wtdiscord_net_Network 1 */
/* MAP_WtTOGGLE_FN(map_toggle__wtdiscord_net_Network){ */
/*   MAP_WtTOGGLE_PROPAGATE; */
/* } */

/* #define map_toggle_maxtoggles__wtintersect_net_Network 1 */
/* MAP_WtTOGGLE_FN(map_toggle__wtintersect_net_Network){ */
/*   ModelTerm *mtp = auxnet->mtp; */
/*   int *ref_el = IINPUT_PARAM; */
/*   MAP_WtTOGGLE_PROPAGATE_IF(iEdgeListSearch(tail, head, ref_el)); */
/* } */

/* #define map_toggle_maxtoggles__wtunion_net_Network 1 */
/* MAP_WtTOGGLE_FN(map_toggle__wtunion_net_Network){ */
/*   ModelTerm *mtp = auxnet->mtp; */
/*   int *ref_el = IINPUT_PARAM; */
/*   MAP_WtTOGGLE_PROPAGATE_IF(!iEdgeListSearch(tail, head, ref_el)); */
/* } */

/* #define map_toggle_maxtoggles__wtblockdiag_net 1 */
/* MAP_WtTOGGLE_FN(map_toggle__wtblockdiag_net){ */
/*   ModelTerm *mtp = auxnet->mtp; */
/*   int *b = IINPUT_PARAM-1; // tail and head are indexed from 1. */
/*   MAP_WtTOGGLE_PROPAGATE_IF(b[tail]==b[head]); */
/* } */

#define __wtundir_net_totoggle				\
  Rboolean totoggle;					\
  double ht = WtGETWT(head,tail);                       \
  double ostate = WtGETWT(head,tail, auxnet->onwp);     \
  double w;                                             \
  switch(rule){						\
  case 1: /* max */					\
    w = MAX(ht, weight);                                \
    totoggle = w != ostate;                             \
    break;						\
  case 2: /* min */					\
    w = MIN(ht, weight);                                \
    totoggle = w != ostate;                             \
    break;						\
  case 3: /* upper */					\
    w = weight;                                         \
    totoggle = tail<=head && w != ostate;               \
    break;						\
  case 4: /* lower */					\
    w = weight;                                         \
    totoggle = tail>=head && w != ostate;               \
    break;						\
  default: /* never reached, but avoids a warning */	\
    w = 0;                                              \
    totoggle = FALSE;					\
  }

#define map_toggle_maxtoggles__wtundir_net 1
MAP_WtTOGGLE_FN(map_toggle__wtundir_net){
  WtModelTerm *mtp = auxnet->mtp;
  unsigned int rule = IINPUT_PARAM[0];
  WtNetwork *nwp = auxnet->inwp;
  __wtundir_net_totoggle;
  if(totoggle){
    *tails = MIN(tail,head);
    *heads = MAX(tail,head);
    *weights = w;
    return 1;
  }else{
    return 0;
  }
}

/* #define __wtfilter_formula_net_totoggle(t, h, w, nwp, m)                 \ */
/*   Rboolean totoggle;                                                    \ */
/*   {                                                                     \ */
/*     WtChangeStats1((t), (h), (w), (nwp), (m), edgestate);               \ */
/*     double change = edgestate ? -*(m->workspace) : *(m->workspace);     \ */
/*     switch(op){                                                         \ */
/*     case 1: totoggle = change == 0; break;                              \ */
/*     case 2: totoggle = change == INPUT_PARAM[0]; break;                 \ */
/*     case 3: totoggle = change != INPUT_PARAM[0]; break;                 \ */
/*     case 4: totoggle = change > INPUT_PARAM[0]; break;                  \ */
/*     case 5: totoggle = change < INPUT_PARAM[0]; break;                  \ */
/*     case 6: totoggle = change >= INPUT_PARAM[0]; break;                 \ */
/*     case 7: totoggle = change <= INPUT_PARAM[0]; break;                 \ */
/*     default: totoggle = change != 0;                                    \ */
/*     }                                                                   \ */
/*   } */

/* #define map_toggle_maxtoggles__wtfilter_formula_net 1 */
/* MAP_WtTOGGLE_FN(map_toggle__wtfilter_formula_net){ */
/*   WtModelTerm *mtp = auxnet->mtp; */
/*   unsigned int op = IINPUT_PARAM[0]; */
/*   GET_STORAGE(WtModel, m); */
/*   __wtfilter_formula_net_totoggle(tail, head, weight, auxnet->inwp, m); */
/*   MAP_WtTOGGLE_PROPAGATE_IF(totoggle); */
/* } */

#define map_toggle_maxtoggles__wtsubgraph_net 1
MAP_WtTOGGLE_FN(map_toggle__wtsubgraph_net){
  WtNetwork *nwp = auxnet->inwp;
  WtModelTerm *mtp = auxnet->mtp;
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
    *weights = weight;
    return 1;
  }else{
    return 0;
  }
}

#define __wttransformed_net_totoggle                    \
  Rboolean totoggle;					\
  double ostate = WtGETWT(head,tail, auxnet->onwp);     \
  double w;                                             \
  switch(op){						\
  case 1: /* sqrt */					\
    w = sqrt(weight);                                   \
    totoggle = w != ostate;                             \
    break;						\
  default: /* never reached, but avoids a warning */	\
    w = 0;                                              \
    totoggle = FALSE;					\
  }

#define map_toggle_maxtoggles__wttransformed_net 1
MAP_WtTOGGLE_FN(map_toggle__wttransformed_net){
  WtModelTerm *mtp = auxnet->mtp;
  unsigned int op = IINPUT_PARAM[0];
  __wttransformed_net_totoggle;
  if(totoggle){
    *tails = tail;
    *heads = head;
    *weights = w;
    return 1;
  }else{
    return 0;
  }
}

#endif
