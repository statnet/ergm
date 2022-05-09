/*  File src/changestats_internal.c in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2022 Statnet Commons
 */
#include "ergm_edgetree.h"
#include "ergm_changestat.h"

/***************************************************************
 changestats internal
***************************************************************/

/*****************
 changestat: d_b1degree_edgecov
*****************/

CHANGESTAT_FN(d_b1degree_edgecov) { 
  int i, j, k, echange, n1, /* n2, */ edgecovval, min, max, mid, nedges;
  Vertex tail1, head1, tail2, head2, b1deg, d;

  n1 = BIPARTITE;
  /* n2 = N_NODES - BIPARTITE; */
  nedges = (N_INPUT_PARAMS - N_CHANGE_STATS)/2;
  for (i=0; i < N_CHANGE_STATS; i++) 
    CHANGE_STAT[i] = 0.0;  
  FOR_EACH_TOGGLE(i) {
    b1deg = 0;
    tail1 = TAIL(i);
    head1 = HEAD(i);
    echange = IS_OUTEDGE(tail1, head1) ? -1 : 1;
    edgecovval = 0;
    min=N_CHANGE_STATS;          /*  Begin: Determine whether this tie exists in edgecov*/
    max=nedges+min-1;
    while (max >= min) {
      mid = (min + max)/2;
      tail2 = INPUT_PARAM[mid+N_CHANGE_STATS-1];
      head2 = INPUT_PARAM[mid+nedges+N_CHANGE_STATS-1];
      if (tail1<tail2 || ((tail1==tail2)&&(head1<head2))) { /* Move search window down */
        max = mid-1;
      }
      else if (tail1>tail2 || ((tail1==tail2)&&(head1>head2))) { /* Move search window up */
        min = mid+1;
      }
      else {
        edgecovval = 1;
        break;
      }
    }                            /*  End: Determine whether this tie exists in edgecov*/
    if (edgecovval == 1) {       /*  Begin: determine current edgecov-degree of node of interest*/
      for (k=n1+1; k <= N_NODES; k++) 
        if (IS_OUTEDGE(tail1,k)) {
          min=N_CHANGE_STATS;
          max=nedges+min-1;
          while (max >= min) {
            mid = (min + max)/2;
            tail2 = INPUT_PARAM[mid+N_CHANGE_STATS-1];
            head2 = INPUT_PARAM[mid+nedges+N_CHANGE_STATS-1];
            if (tail1<tail2 || ((tail1==tail2)&&(k<head2))) { /* Move search window down */
              max = mid-1;
            }
            else if (tail1>tail2 || ((tail1==tail2)&&(k>head2))) { /* Move search window up */
              min = mid+1;
            }
            else {
              b1deg ++;
              break;
            }
          }
        }                      
    }                            /*  End: determine current edgecov-degree of node of interest*/
    for(j = 0; j < N_CHANGE_STATS; j++) {
      d = (Vertex)(INPUT_PARAM[j]);
      CHANGE_STAT[j] += (b1deg + echange*edgecovval == d) - (b1deg == d);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_b2degree_edgecov
*****************/

CHANGESTAT_FN(d_b2degree_edgecov) { 
  int i, j, k, echange, n1, /* n2, */ edgecovval, min, max, mid, nedges;
  Vertex tail1, head1, tail2, head2, b2deg, d;

  n1 = BIPARTITE;
  /* n2 = N_NODES - BIPARTITE; */
  nedges = (N_INPUT_PARAMS - N_CHANGE_STATS)/2;  
  for (i=0; i < N_CHANGE_STATS; i++)
    CHANGE_STAT[i] = 0.0;  
  FOR_EACH_TOGGLE(i) {
    b2deg = 0;
    tail1 = TAIL(i);
    head1 = HEAD(i);
    echange = IS_OUTEDGE(tail1, head1) ? -1 : 1;
    edgecovval = 0;
    min=N_CHANGE_STATS;          /*  Begin: Determine whether this tie exists in edgecov*/
    max=nedges+min-1;
    while (max >= min) {
      mid = (min + max)/2;
      tail2 = INPUT_PARAM[mid+N_CHANGE_STATS-1];
      head2 = INPUT_PARAM[mid+nedges+N_CHANGE_STATS-1];
      if (tail1<tail2 || ((tail1==tail2)&&(head1<head2))) { /* Move search window down */
        max = mid-1;
      }
      else if (tail1>tail2 || ((tail1==tail2)&&(head1>head2))) { /* Move search window up */
        min = mid+1;
      }
      else {
        edgecovval = 1;
        break;                  
      }
    }                            /*  End: Determine whether this tie exists in edgecov*/
    if (edgecovval == 1) {       /*  Begin: determine current edgecov-degree of node of interest*/
      for (k=1; k <= n1; k++) 
        if (IS_OUTEDGE(k,head1)) {
          min=N_CHANGE_STATS;
          max=nedges+min-1;
          while (max >= min) {
            mid = (min + max)/2;
            tail2 = INPUT_PARAM[mid+N_CHANGE_STATS-1];
            head2 = INPUT_PARAM[mid+nedges+N_CHANGE_STATS-1];
            if (k<tail2 || ((k==tail2)&&(head1<head2))) { /* Move search window down */
              max = mid-1;
            }
            else if (k>tail2 || ((k==tail2)&&(head1>head2))) { /* Move search window up */
              min = mid+1;
            }
            else {
              b2deg ++;
              break;
            }
          }
        }                      
    }                            /*  End: determine current edgecov-degree of node of interest*/
    for(j = 0; j < N_CHANGE_STATS; j++) {
      d = (Vertex)(INPUT_PARAM[j]);
      CHANGE_STAT[j] += (b2deg + echange*edgecovval == d) - (b2deg == d);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}


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


/*****************
 changestat: d_b1mindegree_edgecov
*****************/
CHANGESTAT_FN(d_b1mindegree_edgecov) { 
  int i, j, k, echange, n1, /* n2, */ edgecovval, min, max, mid, nedges;
  Vertex tail1, head1, tail2, head2, b1deg, d;

  n1 = BIPARTITE;
  /* n2 = N_NODES - BIPARTITE; */
  nedges = (N_INPUT_PARAMS - N_CHANGE_STATS)/2;
  for (i=0; i < N_CHANGE_STATS; i++) 
    CHANGE_STAT[i] = 0.0;  
  FOR_EACH_TOGGLE(i) {
    b1deg = 0;
    tail1 = TAIL(i);
    head1 = HEAD(i);
    echange = IS_OUTEDGE(tail1, head1) ? -1 : 1;
    edgecovval = 0;
    min=N_CHANGE_STATS;          /*  Begin: Determine whether this tie exists in edgecov*/
    max=nedges+min-1;
    while (max >= min) {
      mid = (min + max)/2;
      tail2 = INPUT_PARAM[mid+N_CHANGE_STATS-1];
      head2 = INPUT_PARAM[mid+nedges+N_CHANGE_STATS-1];
      if (tail1<tail2 || ((tail1==tail2)&&(head1<head2))) { /* Move search window down */
        max = mid-1;
      }
      else if (tail1>tail2 || ((tail1==tail2)&&(head1>head2))) { /* Move search window up */
        min = mid+1;
      }
      else {
        edgecovval = 1;
        break;
      }
    }                            /*  End: Determine whether this tie exists in edgecov*/
    if (edgecovval == 1) {       /*  Begin: determine current edgecov-degree of node of interest*/
      for (k=n1+1; k <= N_NODES; k++) 
        if (IS_OUTEDGE(tail1,k)) {
          min=N_CHANGE_STATS;
          max=nedges+min-1;
          while (max >= min) {
            mid = (min + max)/2;
            tail2 = INPUT_PARAM[mid+N_CHANGE_STATS-1];
            head2 = INPUT_PARAM[mid+nedges+N_CHANGE_STATS-1];
            if (tail1<tail2 || ((tail1==tail2)&&(k<head2))) { /* Move search window down */
              max = mid-1;
            }
            else if (tail1>tail2 || ((tail1==tail2)&&(k>head2))) { /* Move search window up */
              min = mid+1;
            }
            else {
              b1deg ++;
              break;
            }
          }
        }                      
    }                            /*  End: determine current edgecov-degree of node of interest*/
    for(j = 0; j < N_CHANGE_STATS; j++) {
      d = (Vertex)(INPUT_PARAM[j]);
      CHANGE_STAT[j] += (b1deg + echange*edgecovval >= d) - (b1deg >= d);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}


/*****************
 changestat: d_b2mindegree_edgecov
*****************/

CHANGESTAT_FN(d_b2mindegree_edgecov) { 
  int i, j, k, echange, n1, /* n2, */ edgecovval, min, max, mid, nedges;
  Vertex tail1, head1, tail2, head2, b2deg, d;

  n1 = BIPARTITE;
  /* n2 = N_NODES - BIPARTITE; */
  nedges = (N_INPUT_PARAMS - N_CHANGE_STATS)/2;  
  for (i=0; i < N_CHANGE_STATS; i++)
    CHANGE_STAT[i] = 0.0;  
  FOR_EACH_TOGGLE(i) {
    b2deg = 0;
    tail1 = TAIL(i);
    head1 = HEAD(i);
    echange = IS_OUTEDGE(tail1, head1) ? -1 : 1;
    edgecovval = 0;
    min=N_CHANGE_STATS;          /*  Begin: Determine whether this tie exists in edgecov*/
    max=nedges+min-1;
    while (max >= min) {
      mid = (min + max)/2;
      tail2 = INPUT_PARAM[mid+N_CHANGE_STATS-1];
      head2 = INPUT_PARAM[mid+nedges+N_CHANGE_STATS-1];
      if (tail1<tail2 || ((tail1==tail2)&&(head1<head2))) { /* Move search window down */
        max = mid-1;
      }
      else if (tail1>tail2 || ((tail1==tail2)&&(head1>head2))) { /* Move search window up */
        min = mid+1;
      }
      else {
        edgecovval = 1;
        break;                  
      }
    }                            /*  End: Determine whether this tie exists in edgecov*/
    if (edgecovval == 1) {       /*  Begin: determine current edgecov-degree of node of interest*/
      for (k=1; k <= n1; k++) {
        if (IS_OUTEDGE(k,head1)) {
          min=N_CHANGE_STATS;
          max=nedges+min-1;
          while (max >= min) {
            mid = (min + max)/2;
            tail2 = INPUT_PARAM[mid+N_CHANGE_STATS-1];
            head2 = INPUT_PARAM[mid+nedges+N_CHANGE_STATS-1];
            if (k<tail2 || ((k==tail2)&&(head1<head2))) { /* Move search window down */
              max = mid-1;
            }
            else if (k>tail2 || ((k==tail2)&&(head1>head2))) { /* Move search window up */
              min = mid+1;
            }
            else {
              b2deg ++;
              break;
            }
          }
        }
      }
    }          /*  End: determine current edgecov-degree of node of interest*/
    for(j = 0; j < N_CHANGE_STATS; j++) {
      d = (Vertex)(INPUT_PARAM[j]);
      CHANGE_STAT[j] += (b2deg + echange*edgecovval >= d) - (b2deg >= d);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}


