#include "wtchangestats.h"

/********************  changestats:   C    ***********/

/*****************
 stat: cyclicweights(_max)
*****************/

D_CHANGESTAT_FN(d_cyclicweights_max) = d_from_s;

S_CHANGESTAT_FN(s_cyclicweights_max) { 
  Edge e1, e2;
  Vertex h, t, change, node3;
  
  change=0;
  for (h=1; h <= N_NODES; h++) {
    STEP_THROUGH_OUTEDGES(h, e1, t) {
      double best_path = 0;
      STEP_THROUGH_OUTEDGES(t, e2, node3) { 
	best_path = fmax(best_path, fmin(GETWT(node3,h),GETWT(t,node3)));
      }
      CHANGE_STAT[0] += fmin(best_path, GETWT(h,t));
    }
  }
}

/*****************
 stat: cyclicweights(_sum)
*****************/

D_CHANGESTAT_FN(d_cyclicweights_sum) = d_from_s;

S_CHANGESTAT_FN(s_cyclicweights_sum) { 
  Edge e1, e2;
  Vertex h, t, change, node3;
  
  change=0;
  for (h=1; h <= N_NODES; h++) {
    STEP_THROUGH_OUTEDGES(h, e1, t) {
      double path_strength = 0;
      STEP_THROUGH_OUTEDGES(t, e2, node3) { 
	path_strength += fmin(GETWT(node3,h),GETWT(t,node3));
      }
      CHANGE_STAT[0] += fmin(path_strength, GETWT(h,t));
    }
  }
}

/********************  changestats:   N    ***********/

/*****************
 stat: mutualweights
*****************/
WtD_CHANGESTAT_FN(d_mutualweights_min){
  double change, OLDWT, thweight;
  int i;
  
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i){
    GETOLDWT(i);
    thweight = GETWT(tails[i],heads[i]);
    CHANGE_STAT[0] += fmin(weights[i],thweight) - fmin(OLDWT,thweight);
    SETWT_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_SETWTS(i);
}

/*****************
 stat: mutualweights
*****************/
WtD_CHANGESTAT_FN(d_mutualweights_nabsdiff){
  double change, OLDWT, thweight;
  int i;
  
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i){
    GETOLDWT(i);
    thweight = GETWT(tails[i],heads[i]);
    CHANGE_STAT[0] -= fabs(weights[i]-thweight) - fabs(OLDWT-thweight);
    SETWT_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_SETWTS(i);
}


/********************  changestats:   N    ***********/

/*****************
 stat: nonzero
*****************/
WtD_CHANGESTAT_FN(d_nonzero) {
  double OLDWT;
  int i;
  
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i){
    GETOLDWT(i);
    CHANGE_STAT[0] += (weights[i]!=0) - (OLDWT!=0);
    SETWT_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_SETWTS(i);
}



/********************  changestats:   S    ***********/

/*****************
 stat: sum
*****************/
WtD_CHANGESTAT_FN(d_sum) {
  double OLDWT;
  int i;
  
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i){
    GETOLDWT(i);
    CHANGE_STAT[0] += weights[i]-OLDWT;
    SETWT_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_SETWTS(i);
}

/********************  changestats:   T    ***********/

/*****************
 stat: transitiveweights(_max)
*****************/

D_CHANGESTAT_FN(d_transitiveweights_max) = d_from_s;

S_CHANGESTAT_FN(s_transitiveweights_max) { 
  Edge e1, e2;
  Vertex h, t, change, node3;
  
  change=0;
  for (h=1; h <= N_NODES; h++) {
    STEP_THROUGH_OUTEDGES(h, e1, t) {
      double best_path = 0;
      STEP_THROUGH_INEDGES(t, e2, node3) { 
	best_path = fmax(best_path, fmin(GETWT(h,node3),GETWT(node3,t)));
      }
      CHANGE_STAT[0] += fmin(best_path, GETWT(h,t));
    }
  }
}

/*****************
 stat: transitiveweights(_sum)
*****************/

D_CHANGESTAT_FN(d_transitiveweights_sum) = d_from_s;

S_CHANGESTAT_FN(s_transitiveweights_sum) { 
  Edge e1, e2;
  Vertex h, t, change, node3;
  
  change=0;
  for (h=1; h <= N_NODES; h++) {
    STEP_THROUGH_OUTEDGES(h, e1, t) {
      double path_strength = 0;
      STEP_THROUGH_INEDGES(t, e2, node3) { 
	path_strength += fmin(GETWT(h,node3),GETWT(node3,t));
      }
      CHANGE_STAT[0] += fmin(path_strength, GETWT(h,t));
    }
  }
}
