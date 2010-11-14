#include "wtchangestats.h"

/********************  changestats:   C    ***********/

/*****************
 stat: cyclicweights(_max)
*****************/

WtD_FROM_S_FN(d_cyclicweights_max)

WtS_CHANGESTAT_FN(s_cyclicweights_max) { 
  Edge e1, e2;
  Vertex h, t, node3;
  
  CHANGE_STAT[0]=0;
  for (h=1; h <= N_NODES; h++) {
    EXEC_THROUGH_FOUTEDGES(h, e1, t, {
      double best_path = 0;
      EXEC_THROUGH_OUTEDGES(t, e2, node3, { 
	best_path = fmax(best_path, fmin(GETWT(node3,h),GETWT(t,node3)));
	})
      CHANGE_STAT[0] += fmin(best_path, GETWT(h,t));
      })
  }
}

/*****************
 stat: cyclicweights(_sum)
*****************/

WtD_FROM_S_FN(d_cyclicweights_sum)

WtS_CHANGESTAT_FN(s_cyclicweights_sum) { 
  Edge e1, e2;
  Vertex h, t, node3;
  
  CHANGE_STAT[0]=0;
  for (h=1; h <= N_NODES; h++) {
    EXEC_THROUGH_FOUTEDGES(h, e1, t, {
      double path_strength = 0;
      EXEC_THROUGH_OUTEDGES(t, e2, node3, { 
	path_strength += fmin(GETWT(node3,h),GETWT(t,node3));
	})
      CHANGE_STAT[0] += fmin(path_strength, GETWT(h,t));
      })
  }
}

/********************  changestats:   M    ***********/

/*****************
 stat: mutual (minimum)
*****************/
WtD_CHANGESTAT_FN(d_mutual_wt_min){
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
 stat: mutual (-|yij-yji|)
*****************/
WtD_CHANGESTAT_FN(d_mutual_wt_nabsdiff){
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
 stat: nodefactor
*****************/
WtD_CHANGESTAT_FN(d_nodefactor_wt){ 
  double s, factorval, OLDWT;
  Vertex h, t;
  int i, j, hattr, tattr;
  
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    h = heads[i];
    t = tails[i];
    GETOLDWT(i);
    s = weights[i] - OLDWT;
    hattr = INPUT_ATTRIB[h-1];
    tattr = INPUT_ATTRIB[t-1];
    for (j=0; j < N_CHANGE_STATS; j++) {
      factorval = INPUT_PARAM[j];
      if (hattr == factorval) CHANGE_STAT[j] += s;
      if (tattr == factorval) CHANGE_STAT[j] += s;
    }
    SETWT_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_SETWTS(i);
}

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
 stat: nsumlogfactorial
 Turns a Poisson-reference or geometric-reference ERGM into a Conway-Maxwell-Poisson distribution
*****************/
WtD_CHANGESTAT_FN(d_nsumlogfactorial) {
  double OLDWT;
  int i;
  
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i){
    GETOLDWT(i);
    CHANGE_STAT[0] -= lgamma1p(weights[i])-lgamma1p(OLDWT);
    SETWT_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_SETWTS(i);
}

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

WtD_FROM_S_FN(d_transitiveweights_max)

WtS_CHANGESTAT_FN(s_transitiveweights_max) { 
  Edge e1, e2;
  Vertex h, t, node3;
  
  CHANGE_STAT[0]=0;
  for (h=1; h <= N_NODES; h++) {
    EXEC_THROUGH_FOUTEDGES(h, e1, t, {
      double best_path = 0;
      EXEC_THROUGH_INEDGES(t, e2, node3, { 
	best_path = fmax(best_path, fmin(GETWT(h,node3),GETWT(node3,t)));
	})
      CHANGE_STAT[0] += fmin(best_path, GETWT(h,t));
      })
  }
}

/*****************
 stat: transitiveweights(_sum)
*****************/

WtD_FROM_S_FN(d_transitiveweights_sum)

WtS_CHANGESTAT_FN(s_transitiveweights_sum) { 
  Edge e1, e2;
  Vertex h, t, node3;
  
  CHANGE_STAT[0]=0;
  for (h=1; h <= N_NODES; h++) {
    EXEC_THROUGH_FOUTEDGES(h, e1, t, {
      double path_strength = 0;
      EXEC_THROUGH_INEDGES(t, e2, node3, { 
	path_strength += fmin(GETWT(h,node3),GETWT(node3,t));
	})
      CHANGE_STAT[0] += fmin(path_strength, GETWT(h,t));
      })
  }
}
