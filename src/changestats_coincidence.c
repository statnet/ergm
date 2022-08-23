/*  File src/changestats_coincidence.c in package ergm, part of the
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

/*****************
 changestat: d_coincidence
*****************/
C_CHANGESTAT_FN(c_coincidence) { 
  /*The inputparams are assumed to be set up as follows:
    The first 2*nstats values are in pairs:  (degree, attrvalue)
    The values following the first 2*nstats values are the nodal attributes.*/
  int echange;
  Vertex b1, b2, index;
  Vertex film1, film2, act;
  Vertex nb, nb1, nb2;
  Edge e;
  nb1 = BIPARTITE;
  nb2 = N_NODES - BIPARTITE;
  nb  = 2*nb2-1;

  /* *** don't forget act -> film1 */    
    echange = IS_OUTEDGE(act=tail, film1=head) ? -1 : 1;
    b1 = film1-nb1;
    STEP_THROUGH_OUTEDGES(act, e, film2) {
     if(film2 != film1){
      b2 = film2-nb1;
      if(film2 < film1){
       index = (b2*(nb-b2))/2+b1-nb2;
      }else{
       index = (b1*(nb-b1))/2+b2-nb2;
      }
     if(INPUT_PARAM[index-1]>0.0) CHANGE_STAT[(int)(INPUT_PARAM[index-1])-1] += echange;
     }
    }
}
