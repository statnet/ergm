#include "changestats_transitiveties.h"

/* d_transitiveties and s_transitiveties have been modified and added
   to changestats.c.  The versions here are obsolete (e.g., they still
   assume heads->tails) */
//D_CHANGESTAT_FN(d_transitiveties) { 
//  int i;
//  double current;
//
//  (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */
//  current = mtp->dstats[0];
//  FOR_EACH_TOGGLE(i) { TOGGLE(heads[i],tails[i]); }
//  (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */
//  mtp->dstats[0] -= current;
//  FOR_EACH_TOGGLE(i) { TOGGLE(heads[i],tails[i]); }
//}
//
//
///*****************
// globalstat: s_transitiveties
//*****************/
//S_CHANGESTAT_FN(s_transitiveties) { 
//  Edge e1, e2;
//  Vertex h, t, change, node3;
//  double hattr;
//  int hnottrans;
//  
//  change=0;
//  if(N_INPUT_PARAMS > 0){ /* match on attributes */
//    for (h=1; h <= N_NODES; h++) { 
//      hattr = INPUT_ATTRIB[h-1];
//      STEP_THROUGH_OUTEDGES(h, e1, t) {
//       if(hattr == INPUT_ATTRIB[t-1]) {
//        hnottrans=1;
//          STEP_THROUGH_INEDGES(t, e2, node3) { 
//          if(hnottrans && IS_INEDGE(node3, h) && (hattr == INPUT_ATTRIB[node3-1])){ /* h -> t base forms transitive */
//	   hnottrans=0;
//           change++;
//	  }
//        }
//       }
//      }
//    }
//  }else{
//  for (h=1; h <= N_NODES; h++) { 
//    STEP_THROUGH_OUTEDGES(h, e1, t) {
//      hnottrans=1;
//      STEP_THROUGH_INEDGES(t, e2, node3) { 
//        if(hnottrans && IS_INEDGE(node3, h)){ /* h -> t base forms transitive */
//	 hnottrans=0;
//         change++;
//	}
//      }
//    }
//  }
//  }
//  CHANGE_STAT[0] = change;
//}
D_CHANGESTAT_FN(d_transitiveties2) { 
  int i;
  double current;

  (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */
  current = mtp->dstats[0];
  FOR_EACH_TOGGLE(i) { TOGGLE(heads[i],tails[i]); }
  (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */
  mtp->dstats[0] -= current;
  FOR_EACH_TOGGLE(i) { TOGGLE(heads[i],tails[i]); }
}
/*****************
 globalstat: s_transitiveties2
*****************/
S_CHANGESTAT_FN(s_transitiveties2) { 
  Edge e1, e2;
  Vertex h, t, change, node3;
  double hattr;
  int hnottrans;
  
  change=0;
  if(N_INPUT_PARAMS > 0){ /* match on attributes */
    for (h=1; h <= N_NODES; h++) { 
      hattr = INPUT_ATTRIB[h-1];
      STEP_THROUGH_OUTEDGES(h, e1, t) {
       if(hattr == INPUT_ATTRIB[t-1]) {
        hnottrans=0;
        STEP_THROUGH_INEDGES(t, e2, node3) { 
          if(IS_INEDGE(node3, h) && (hattr == INPUT_ATTRIB[node3-1])){ /* h -> t base forms transitive */
	   hnottrans++;
	  }
        }
        if(hnottrans>2){hnottrans=2;}
        change += hnottrans;
       }
      }
    }
  }else{
  for (h=1; h <= N_NODES; h++) { 
    STEP_THROUGH_OUTEDGES(h, e1, t) {
      hnottrans=0;
      STEP_THROUGH_INEDGES(t, e2, node3) { 
        if(IS_INEDGE(node3, h)){ /* h -> t base forms transitive */
	 hnottrans++;
	}
      }
      if(hnottrans>2){hnottrans=2;}
      change += hnottrans;
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
  FOR_EACH_TOGGLE(i) { TOGGLE(heads[i],tails[i]); }
  (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */
  mtp->dstats[0] -= current;
  FOR_EACH_TOGGLE(i) { TOGGLE(heads[i],tails[i]); }
}
/*****************
 globalstat: s_cyclicalties
*****************/
S_CHANGESTAT_FN(s_cyclicalties) { 
  Edge e1, e2;
  Vertex h, t, change, node3;
  double hattr;
  int hnottrans;
  
  change=0;
  if(N_INPUT_PARAMS > 0){ /* match on attributes */
    for (h=1; h <= N_NODES; h++) { 
      hattr = INPUT_ATTRIB[h-1];
      STEP_THROUGH_OUTEDGES(h, e1, t) {
       if(hattr == INPUT_ATTRIB[t-1]) {
        hnottrans=0;
        STEP_THROUGH_INEDGES(h, e2, node3) { 
          if(IS_INEDGE(node3, t) && (hattr == INPUT_ATTRIB[node3-1])){ /* h -> t base forms cyclical */
	   hnottrans++;
	  }
        }
        if(hnottrans>1){hnottrans=1;}
        change += hnottrans;
       }
      }
    }
  }else{
  for (h=1; h <= N_NODES; h++) { 
    STEP_THROUGH_OUTEDGES(h, e1, t) {
      hnottrans=0;
      STEP_THROUGH_INEDGES(h, e2, node3) { 
        if(IS_INEDGE(node3, t)){ /* h -> t base forms cyclical */
	 hnottrans++;
	}
      }
      if(hnottrans>1){hnottrans=1;}
      change += hnottrans;
    }
  }
  }
  CHANGE_STAT[0] = change;
}
D_CHANGESTAT_FN(d_cyclicalties) { 
  int i;
  double current;

  (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */
  current = mtp->dstats[0];
  FOR_EACH_TOGGLE(i) { TOGGLE(heads[i],tails[i]); }
  (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */
  mtp->dstats[0] -= current;
  FOR_EACH_TOGGLE(i) { TOGGLE(heads[i],tails[i]); }
}
/*****************
 globalstat: s_cyclicalties2
*****************/
S_CHANGESTAT_FN(s_cyclicalties2) { 
  Edge e1, e2;
  Vertex h, t, change, node3;
  double hattr;
  int hnottrans;
  
  change=0;
  if(N_INPUT_PARAMS > 0){ /* match on attributes */
    for (h=1; h <= N_NODES; h++) { 
      hattr = INPUT_ATTRIB[h-1];
      STEP_THROUGH_OUTEDGES(h, e1, t) {
       if(hattr == INPUT_ATTRIB[t-1]) {
        hnottrans=0;
        STEP_THROUGH_INEDGES(h, e2, node3) { 
          if(IS_INEDGE(node3, t) && (hattr == INPUT_ATTRIB[node3-1])){ /* h -> t base forms cyclical */
	   hnottrans++;
	  }
        }
        if(hnottrans>2){hnottrans=2;}
        change += hnottrans;
       }
      }
    }
  }else{
  for (h=1; h <= N_NODES; h++) { 
    STEP_THROUGH_OUTEDGES(h, e1, t) {
      hnottrans=0;
      STEP_THROUGH_INEDGES(h, e2, node3) { 
        if(IS_INEDGE(node3, t)){ /* h -> t base forms cyclical */
	 hnottrans++;
	}
      }
      if(hnottrans>2){hnottrans=2;}
      change += hnottrans;
    }
  }
  }
  CHANGE_STAT[0] = change;
}
///*****************
// changestat: d_transitiveties
//*****************/
//D_CHANGESTAT_FN(d_transitiveties) { 
//  Edge e1, e2, e3;
//  Vertex h, t, change, node3, node4, node5;
//  int i, j, htrans, ntrans, edgeval;
//  double hattr, edgemult;
//  
////Rprintf("call\n");
////  for (node3=1; node3 <= N_NODES; node3++) { 
////  for (node4=1; node4 <= N_NODES; node4++) { 
////   if(IS_OUTEDGE(node3, node4)){Rprintf("h %d t %d\n",node3, node4);}
////  }}
//  ZERO_ALL_CHANGESTATS(i);
//  FOR_EACH_TOGGLE(i) {
//    h = heads[i];
//    t = tails[i];
//    edgeval = IS_OUTEDGE(h, t);
//    htrans=0;
//    for (node3=1; node3 <= N_NODES; node3++) { 
//     if(IS_OUTEDGE(node3, t)&&IS_OUTEDGE(h, node3) ){htrans++;}
//    }
//Rprintf("h %d t %d edgeval %d trit %d\n",h,t,edgeval,htrans);
//    edgemult = edgeval ? -1.0 : 1.0;
//    change = 0;
//    if(N_INPUT_PARAMS > 0){ /* match on attributes */
//      hattr = INPUT_ATTRIB[h-1];
//      if(hattr == INPUT_ATTRIB[t-1]) {
//        STEP_THROUGH_OUTEDGES(t, e1, node3) { /* step through outedges of tail */
//          if(hattr == INPUT_ATTRIB[node3-1])
//            change += IS_INEDGE(node3, h);
//        }
//        STEP_THROUGH_INEDGES(t, e1, node3) { /* step through inedges of tail */
//          if(hattr == INPUT_ATTRIB[node3-1])
//            change += IS_OUTEDGE(node3, h) + IS_INEDGE(node3, h);
//        }
//        if(N_CHANGE_STATS > 1) { /* diff = TRUE; matches must be tabled */
//          for (j=0; j<N_CHANGE_STATS; j++){
//            if (hattr == INPUT_PARAM[j])
//              CHANGE_STAT[j] += edgemult * change;
//          }
//        } else { /* diff = FALSE; all matches equivalent */
//              CHANGE_STAT[0] += edgemult * change;          
//        }
//      }
//    }else{ /* no attribute matching */
//      change=0;
//      ntrans=0;
//      htrans=0;
//      STEP_THROUGH_INEDGES(t, e1, node3) { 
//        if(IS_INEDGE(node3, h)){ /* h -> t base forms transitive, so check */
//         htrans=1;
//	}
//      }
//      if(htrans){
//Rprintf("trans h %d t %d edgeval %d\n",h,t,edgeval);
//       ntrans=0;
//       STEP_THROUGH_OUTEDGES(h, e1, node4) { 
//        if(IS_INEDGE(t, node4)){ ntrans++; }
//       }
//       if(edgeval){
//        if(ntrans==1){change--;}
////	  if(ntrans==1){
////Rprintf("- h %d t %d edgeval %d ntrans %d change %d\n",h,t,edgeval,ntrans,change);
////}
//       }else{
//        if(ntrans==0){change++;}
////	  if(ntrans==0){
////Rprintf("h %d t %d edgeval %d ntrans %d change %d\n",h,t,edgeval,ntrans,change);
////}
//       }
//      }
//      ntrans=0;
//      STEP_THROUGH_INEDGES(t, e1, node3) { 
//        if(IS_OUTEDGE(node3, h)){ /* node3 -> t base forms transitive, so check  */
//         ntrans=0;
//         STEP_THROUGH_OUTEDGES(node3, e2, node4) { 
//           if(IS_OUTEDGE(node4, t)){ ntrans++; }
//	 }
//         if(edgeval){
//          if(ntrans==1){change--;}
////	  if(ntrans==1){
////Rprintf("- node3 %d h %d t %d edgeval %d ntrans %d change %d\n",node3,h,t,edgeval,ntrans,change);
////}
//         }else{
//          if(ntrans==0){change++;}
////	  if(ntrans==0){
////Rprintf("node3 %d h %d t %d edgeval %d ntrans %d change %d\n",node3,h,t,edgeval,ntrans,change);
////}
//         }
//	}
//      }
//      ntrans=0;
//      STEP_THROUGH_OUTEDGES(t, e1, node3) { 
//        if(IS_INEDGE(node3, h)){ /* h -> node3 base forms transitive, so check */
//         ntrans=0;
//         STEP_THROUGH_INEDGES(node3, e2, node4) { 
//           if(IS_OUTEDGE(h, node4)){ ntrans++; }
//	 }
//         if(edgeval){
//          if(ntrans==1){change--;}
////	  if(ntrans==1){
////Rprintf("- node3 %d h %d t %d edgeval %d ntrans %d change %d\n",node3,h,t,edgeval,ntrans,change);
////}
//         }else{
//          if(ntrans==0){change++;}
////	  if(ntrans==0){
////Rprintf("node3 %d h %d t %d edgeval %d ntrans %d change %d\n",node3,h,t,edgeval,ntrans,change);
////}
//         }
//	}
//      }
//      CHANGE_STAT[0] += change;
//      Rprintf("h %d t %d edgeval %d ntrans %d change %d\n",h,t,edgeval,ntrans,change);
//    }
//    TOGGLE_IF_MORE_TO_COME(i);
//  }
////Rprintf("changestat %f change %d\n",CHANGE_STAT[0],change);
//  UNDO_PREVIOUS_TOGGLES(i);
//}
