/*  File src/changestats_internal.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */
#include "ergm_edgetree.h"
#include "ergm_changestat.h"

/*****************
 changestat: d_b1mindegree
*****************/
CHANGESTAT_FN(d_b1mindegree) { 
  int i, j, echange;
  Vertex b1, b1deg, d;

  for (i=0; i < N_CHANGE_STATS; i++) 
    CHANGE_STAT[i] = 0.0;  
  FOR_EACH_TOGGLE(i) {
    b1 = TAIL(i);
    echange = IS_OUTEDGE(b1, HEAD(i)) ? -1 : 1;
    b1deg = OUT_DEG[b1];
    for(j = 0; j < N_CHANGE_STATS; j++) {
      d = (Vertex)(INPUT_PARAM[j]);
      CHANGE_STAT[j] += (b1deg + echange >= d) - (b1deg >= d);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_b2mindegree
*****************/
CHANGESTAT_FN(d_b2mindegree) { 
  /* It is assumed that in this bipartite network, the only edges are
  of the form (b1, b2), where b1 is always strictly less
  than b2.  In other words, the degree of a b1 is equivalent
  to its outdegree and the degree of a b2 is equivalent to its
  indegree.
  */
  int i, j, echange;
  Vertex b1, b2, b2deg, d, *id;

  id=IN_DEG;
  for (i=0; i < N_CHANGE_STATS; i++) 
    CHANGE_STAT[i] = 0.0;  
  for (i=0; i<ntoggles; i++) {      
    echange=IS_OUTEDGE(b1=TAIL(i), b2=HEAD(i)) ? -1 : 1;
    b2deg = id[b2];
    for(j = 0; j < N_CHANGE_STATS; j++) {
      d = (Vertex)(INPUT_PARAM[j]);
      CHANGE_STAT[j] += (b2deg + echange >= d) - (b2deg >= d);
    }
    TOGGLE_IF_MORE_TO_COME(i); 
  }
  UNDO_PREVIOUS_TOGGLES(i);
}
