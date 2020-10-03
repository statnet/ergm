/*  File src/changestats_test.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2020 Statnet Commons
 */
#include "ergm_changestat.h"

/*****************
 changestat: d__edges_times
*****************/
D_CHANGESTAT_FN(d__edges_times) {
  int edgeflag, i;
  Vertex tail, head;
  
  /* *** don't forget tail -> head */    
  CHANGE_STAT[0] = 0.0;
  for (i=0; i < ntoggles; i++)
    {
      edgeflag = IS_OUTEDGE(tail = TAIL(i), head = HEAD(i));
      CHANGE_STAT[0] += edgeflag ? - *INPUT_PARAM : *INPUT_PARAM;
      TOGGLE_IF_MORE_TO_COME(i);
    }
  UNDO_PREVIOUS_TOGGLES(i);
}

S_CHANGESTAT_FN(s__edges_tests) {
  CHANGE_STAT[0] = N_EDGES * *INPUT_PARAM;
}
