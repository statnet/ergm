/*  File src/changestats_indices.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2020 Statnet Commons
 */
#include "changestats_indices.h"

/*****************                       
 changestat: d_indices
*****************/
D_CHANGESTAT_FN(d_indices) { 
  int i;

  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    int t = tails[i];
    int h = heads[i];
    int edgeflag = IS_OUTEDGE(t,h);
    CHANGE_STAT[0] += edgeflag ? -t : t;
    CHANGE_STAT[1] += edgeflag ? -h : h;
    TOGGLE_IF_MORE_TO_COME(i); /* Needed in case of multiple toggles */
  }
  UNDO_PREVIOUS_TOGGLES(i); /* Needed on exit in case of multiple toggles */
}
