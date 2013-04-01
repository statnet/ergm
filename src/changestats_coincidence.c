#include "changestats_coincidence.h"

/*****************
 changestat: d_coincidence
*****************/
D_CHANGESTAT_FN(d_coincidence) { 
  /*The inputparams are assumed to be set up as follows:
    The first 2*nstats values are in pairs:  (degree, attrvalue)
    The values following the first 2*nstats values are the nodal attributes.*/
  int i, echange;
  Vertex b1, b2, index;
  Vertex film1, film2, act;
  Vertex nb, nb1, nb2;
  Edge e;
  nb1 = BIPARTITE;
  nb2 = N_NODES - BIPARTITE;
  nb  = 2*nb2-1;

  /* *** don't forget act -> film1 */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    echange = IS_OUTEDGE(act=TAIL(i), film1=HEAD(i)) ? -1 : 1;
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
    TOGGLE_IF_MORE_TO_COME(i); /* Needed in case of multiple toggles */
  }
  UNDO_PREVIOUS_TOGGLES(i); /* Needed on exit in case of multiple toggles */
}
