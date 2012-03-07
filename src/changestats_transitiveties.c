/*
 *  File ergm/src/changestats_transitiveties.c
 *  Part of the statnet package, http://statnet.org
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) in
 *    http://statnet.org/attribution
 *
 *  Copyright 2012 the statnet development team
 */
#include "changestats_transitiveties.h"

D_CHANGESTAT_FN(d_transitiveties2) { 
  int i;
  double current;

  (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */
  current = mtp->dstats[0];
  FOR_EACH_TOGGLE(i) { TOGGLE(TAIL(i),HEAD(i)); }
  (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */
  mtp->dstats[0] -= current;
  FOR_EACH_TOGGLE(i) { TOGGLE(TAIL(i),HEAD(i)); }
}
/*****************
 globalstat: s_transitiveties2
*****************/
S_CHANGESTAT_FN(s_transitiveties2) { 
  Edge e1, e2;
  Vertex tail, head, change, node3;
  double tailattr;
  int tailnottrans;
  
  change=0;
  if(N_INPUT_PARAMS > 0){ /* match on attributes */
    for (tail=1; tail <= N_NODES; tail++) { 
      tailattr = INPUT_ATTRIB[tail-1];
      STEP_THROUGH_OUTEDGES(tail, e1, head) {
        if(tailattr == INPUT_ATTRIB[head-1]) {
          tailnottrans=0;
          STEP_THROUGH_INEDGES(head, e2, node3) { 
            if(IS_INEDGE(node3, tail) && (tailattr == INPUT_ATTRIB[node3-1])){ /* tail -> head base forms transitive */
              tailnottrans++;
            }
          }
          if(tailnottrans>2){tailnottrans=2;}
          change += tailnottrans;
        }
      }
    }
  }else{
    for (tail=1; tail <= N_NODES; tail++) { 
      STEP_THROUGH_OUTEDGES(tail, e1, head) {
        tailnottrans=0;
        STEP_THROUGH_INEDGES(head, e2, node3) { 
          if(IS_INEDGE(node3, tail)){ /* tail -> head base forms transitive */
            tailnottrans++;
          }
        }
        if(tailnottrans>2){tailnottrans=2;}
        change += tailnottrans;
      }
    }
  }
  CHANGE_STAT[0] = change;
}
D_CHANGESTAT_FN(d_cyclicalties2) { 
  int i;
  double current;

  (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */
  current = mtp->dstats[0];
  FOR_EACH_TOGGLE(i) { TOGGLE(TAIL(i),HEAD(i)); }
  (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */
  mtp->dstats[0] -= current;
  FOR_EACH_TOGGLE(i) { TOGGLE(TAIL(i),HEAD(i)); }
}
///*****************
// globalstat: s_cyclicalties
//*****************/
//S_CHANGESTAT_FN(s_cyclicalties) { 
//  Edge e1, e2;
//  Vertex tail, head, change, node3;
//  double tailattr;
//  int tailnottrans;
//  
//  change=0;
//  if(N_INPUT_PARAMS > 0){ /* match on attributes */
//    for (tail=1; tail <= N_NODES; tail++) { 
//      tailattr = INPUT_ATTRIB[tail-1];
//      STEP_THROUGH_OUTEDGES(tail, e1, head) {
//        if(tailattr == INPUT_ATTRIB[head-1]) {
//          tailnottrans=0;
//          STEP_THROUGH_INEDGES(tail, e2, node3) { 
//            if(IS_INEDGE(node3, head) && (tailattr == INPUT_ATTRIB[node3-1])){ /* tail -> head base forms cyclical */
//              tailnottrans++;
//            }
//          }
//          if(tailnottrans>1){tailnottrans=1;}
//          change += tailnottrans;
//        }
//      }
//    }
//  }else{
//    for (tail=1; tail <= N_NODES; tail++) { 
//      STEP_THROUGH_OUTEDGES(tail, e1, head) {
//        tailnottrans=0;
//        STEP_THROUGH_INEDGES(tail, e2, node3) { 
//          if(IS_INEDGE(node3, head)){ /* tail -> head base forms cyclical */
//            tailnottrans++;
//          }
//        }
//        if(tailnottrans>1){tailnottrans=1;}
//        change += tailnottrans;
//      }
//    }
//  }
//  CHANGE_STAT[0] = change;
//}
//D_CHANGESTAT_FN(d_cyclicalties) { 
//  int i;
//  double current;
//
//  (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */
//  current = mtp->dstats[0];
//  FOR_EACH_TOGGLE(i) { TOGGLE(TAIL(i),HEAD(i)); }
//  (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */
//  mtp->dstats[0] -= current;
//  FOR_EACH_TOGGLE(i) { TOGGLE(TAIL(i),HEAD(i)); }
//}
/*****************
 globalstat: s_cyclicalties2
*****************/
S_CHANGESTAT_FN(s_cyclicalties2) { 
  Edge e1, e2;
  Vertex tail, head, change, node3;
  double tailattr;
  int tailnottrans;
  
  change=0;
  if(N_INPUT_PARAMS > 0){ /* match on attributes */
    for (tail=1; tail <= N_NODES; tail++) { 
      tailattr = INPUT_ATTRIB[tail-1];
      STEP_THROUGH_OUTEDGES(tail, e1, head) {
        if(tailattr == INPUT_ATTRIB[head-1]) {
          tailnottrans=0;
          STEP_THROUGH_INEDGES(tail, e2, node3) { 
            if(IS_INEDGE(node3, head) && (tailattr == INPUT_ATTRIB[node3-1])){ /* tail -> head base forms cyclical */
              tailnottrans++;
            }
          }
          if(tailnottrans>2){tailnottrans=2;}
          change += tailnottrans;
        }
      }
    }
  }else{
    for (tail=1; tail <= N_NODES; tail++) { 
      STEP_THROUGH_OUTEDGES(tail, e1, head) {
        tailnottrans=0;
        STEP_THROUGH_INEDGES(tail, e2, node3) { 
          if(IS_INEDGE(node3, head)){ /* tail -> head base forms cyclical */
            tailnottrans++;
          }
        }
        if(tailnottrans>2){tailnottrans=2;}
        change += tailnottrans;
      }
    }
  }
  CHANGE_STAT[0] = change;
}

