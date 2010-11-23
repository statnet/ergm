#include "wtchangestats.h"

/********************  changestats:   A    ***********/

/*****************                       
 stat: absdiff(_nonzero)
*****************/
WtD_CHANGESTAT_FN(d_absdiff_nonzero){ 
  double p = INPUT_ATTRIB[0];
  int i;

  ZERO_ALL_CHANGESTATS(i);
  EXEC_THROUGH_TOGGLES(i,{
    if(p==1.0){
      CHANGE_STAT[0] += fabs(INPUT_ATTRIB[heads[i]] - INPUT_ATTRIB[tails[i]])*((weights[i]!=0)-(CURWT!=0));
    } else {
      CHANGE_STAT[0] += pow(fabs(INPUT_ATTRIB[heads[i]] - INPUT_ATTRIB[tails[i]]), p)*((weights[i]!=0)-(CURWT!=0));
    }
    })
}

/*****************                       
 stat: absdiff(_sum)
*****************/
WtD_CHANGESTAT_FN(d_absdiff_sum){ 
  double p = INPUT_ATTRIB[0];
  int i;

  ZERO_ALL_CHANGESTATS(i);
  EXEC_THROUGH_TOGGLES(i,{
    if(p==1.0){
      CHANGE_STAT[0] += fabs(INPUT_ATTRIB[heads[i]] - INPUT_ATTRIB[tails[i]])*(weights[i]-CURWT);
    } else {
      CHANGE_STAT[0] += pow(fabs(INPUT_ATTRIB[heads[i]] - INPUT_ATTRIB[tails[i]]), p)*(weights[i]-CURWT);
    }
    })
}

/*****************
 stat: absdiffcat(_nonzero)
*****************/
WtD_CHANGESTAT_FN(d_absdiffcat_nonzero){ 
  double change, absdiff, NAsubstitute, hval, tval;
  Vertex ninputs;
  int i, j;
  
  ninputs = N_INPUT_PARAMS - N_NODES;
  NAsubstitute = INPUT_PARAM[ninputs-1];
  ZERO_ALL_CHANGESTATS(i);
  EXEC_THROUGH_TOGGLES(i,{
      change = (weights[i]!=0)-(CURWT!=0);
      hval = INPUT_ATTRIB[heads[i]-1];
      tval = INPUT_ATTRIB[tails[i]-1];
      if (hval == NAsubstitute ||  tval == NAsubstitute) absdiff = NAsubstitute;
      else absdiff = fabs(hval - tval);
      if (absdiff>0){
	for (j=0; j<N_CHANGE_STATS; j++){
	  CHANGE_STAT[j] += (absdiff==INPUT_PARAM[j]) ? change : 0.0;
	}
      }
    })
}

/*****************
 stat: absdiffcat(_sum)
*****************/
WtD_CHANGESTAT_FN(d_absdiffcat_sum){ 
  double change, absdiff, NAsubstitute, hval, tval;
  Vertex ninputs;
  int i, j;
  
  ninputs = N_INPUT_PARAMS - N_NODES;
  NAsubstitute = INPUT_PARAM[ninputs-1];
  ZERO_ALL_CHANGESTATS(i);
  EXEC_THROUGH_TOGGLES(i,{
      change = weights[i]-CURWT;
      hval = INPUT_ATTRIB[heads[i]-1];
      tval = INPUT_ATTRIB[tails[i]-1];
      if (hval == NAsubstitute ||  tval == NAsubstitute) absdiff = NAsubstitute;
      else absdiff = fabs(hval - tval);
      if (absdiff>0){
	for (j=0; j<N_CHANGE_STATS; j++){
	  CHANGE_STAT[j] += (absdiff==INPUT_PARAM[j]) ? change : 0.0;
	}
      }
    })
}


/*****************
 stat: atleast
*****************/
WtD_CHANGESTAT_FN(d_atleast){
  int i;
  
  CHANGE_STAT[0] = 0.0;
  EXEC_THROUGH_TOGGLES(i,{
      CHANGE_STAT[0] += (weights[i]>=INPUT_ATTRIB[0]) - (CURWT>=INPUT_ATTRIB[0]);
    })
}

/********************  changestats:   C    ***********/

/*****************
 stat: cyclicweights(_max)
*****************/

WtD_FROM_S_FN(d_cyclicweights_max)

WtS_CHANGESTAT_FN(s_cyclicweights_max){ 
  Edge e1, e2;
  Vertex h, t, node3;
  
  CHANGE_STAT[0]=0;
  for (h=1; h <= N_NODES; h++){
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

WtS_CHANGESTAT_FN(s_cyclicweights_sum){ 
  Edge e1, e2;
  Vertex h, t, node3;
  
  CHANGE_STAT[0]=0;
  for (h=1; h <= N_NODES; h++){
    EXEC_THROUGH_FOUTEDGES(h, e1, t, {
      double path_strength = 0;
      EXEC_THROUGH_OUTEDGES(t, e2, node3, { 
	path_strength += fmin(GETWT(node3,h),GETWT(t,node3));
	})
      CHANGE_STAT[0] += fmin(path_strength, GETWT(h,t));
      })
  }
}

/*****************
 stat: cyclicweights(_threshold)
*****************/

WtD_FROM_S_FN(d_cyclicweights_threshold)

WtS_CHANGESTAT_FN(s_cyclicweights_threshold){ 
  Edge e1, e2;
  Vertex h, t, node3;
  
  CHANGE_STAT[0]=0;
  for (h=1; h <= N_NODES; h++){
    EXEC_THROUGH_FOUTEDGES(h, e1, t, {
	if(GETWT(h,t)<=INPUT_ATTRIB[0]) break;
	EXEC_THROUGH_OUTEDGES(t, e2, node3, { 
	    if(GETWT(node3,h)>INPUT_ATTRIB[0] && GETWT(t,node3)>INPUT_ATTRIB[0]){
	      CHANGE_STAT[0]++; 
	      break;
	    }
	  })
	  })
      }
}


/********************  changestats:   G    ***********/

/*****************
 stat: greaterthan
*****************/
WtD_CHANGESTAT_FN(d_greaterthan){
  int i;
  
  CHANGE_STAT[0] = 0.0;
  EXEC_THROUGH_TOGGLES(i,{
      CHANGE_STAT[0] += (weights[i]>INPUT_ATTRIB[0]) - (CURWT>INPUT_ATTRIB[0]);
  })
}

/********************  changestats:   I    ***********/

/*****************
 stat: ininterval
*****************/
WtD_CHANGESTAT_FN(d_ininterval){
  int i;
  
  CHANGE_STAT[0] = 0.0;
  EXEC_THROUGH_TOGGLES(i,{
      CHANGE_STAT[0] += ((INPUT_ATTRIB[2] ? weights[i]>INPUT_ATTRIB[0] : weights[i]>=INPUT_ATTRIB[0]) && (INPUT_ATTRIB[3] ? weights[i]<INPUT_ATTRIB[1] : weights[i]<=INPUT_ATTRIB[1])) - ((INPUT_ATTRIB[2] ? CURWT>INPUT_ATTRIB[0] : CURWT>=INPUT_ATTRIB[0]) && (INPUT_ATTRIB[3] ? CURWT<INPUT_ATTRIB[1] : CURWT<=INPUT_ATTRIB[1]));
    })
}


/********************  changestats:   M    ***********/

/*****************
 stat: mutual (minimum)
*****************/
WtD_CHANGESTAT_FN(d_mutual_wt_min){
  double thweight;
  int i;
  
  ZERO_ALL_CHANGESTATS(i);
  EXEC_THROUGH_TOGGLES(i,{
      thweight = GETWT(tails[i],heads[i]);
      CHANGE_STAT[0] += fmin(weights[i],thweight) - fmin(CURWT,thweight);
    })
}


/*****************
 stat: mutual (-|yij-yji|)
*****************/
WtD_CHANGESTAT_FN(d_mutual_wt_nabsdiff){
  double thweight;
  int i;
  
  ZERO_ALL_CHANGESTATS(i);
  EXEC_THROUGH_TOGGLES(i,{
      thweight = GETWT(tails[i],heads[i]);
    CHANGE_STAT[0] -= fabs(weights[i]-thweight) - fabs(CURWT-thweight);
  })
}

/********************  changestats:   N    ***********/

/*****************
 stat: nodecov (nonzero)
*****************/
WtD_CHANGESTAT_FN(d_nodecov_nonzero){ 
  int i;

  CHANGE_STAT[0] = 0.0;
  EXEC_THROUGH_TOGGLES(i,{
      CHANGE_STAT[0] += (INPUT_ATTRIB[heads[i]-1] + INPUT_ATTRIB[tails[i]-1])*((weights[i]!=0)-(CURWT!=0));
  })
}

/*****************
 stat: d_nodecov
*****************/
WtD_CHANGESTAT_FN(d_nodecov_sum){ 
  int i;

  CHANGE_STAT[0] = 0.0;
  EXEC_THROUGH_TOGGLES(i,{
      CHANGE_STAT[0] += (INPUT_ATTRIB[heads[i]-1] + INPUT_ATTRIB[tails[i]-1])*(weights[i]-CURWT);
  })
}


/*****************
 stat: nodefactor (nonzero)
*****************/
WtD_CHANGESTAT_FN(d_nodefactor_nonzero){ 
  double s, factorval;
  Vertex h, t;
  int i, j, hattr, tattr;
  
  ZERO_ALL_CHANGESTATS(i);
  EXEC_THROUGH_TOGGLES(i, {
    h = heads[i];
    t = tails[i];
        s = (weights[i]!=0) - (CURWT!=0);
    hattr = INPUT_ATTRIB[h-1];
    tattr = INPUT_ATTRIB[t-1];
    for (j=0; j < N_CHANGE_STATS; j++){
      factorval = INPUT_PARAM[j];
      if (hattr == factorval) CHANGE_STAT[j] += s;
      if (tattr == factorval) CHANGE_STAT[j] += s;
    }
  })
}

/*****************
 stat: nodefactor (sum)
*****************/
WtD_CHANGESTAT_FN(d_nodefactor_sum){ 
  double s, factorval;
  Vertex h, t;
  int i, j, hattr, tattr;
  
  ZERO_ALL_CHANGESTATS(i);
  EXEC_THROUGH_TOGGLES(i, {
    h = heads[i];
    t = tails[i];
    s = weights[i] - CURWT;
    hattr = INPUT_ATTRIB[h-1];
    tattr = INPUT_ATTRIB[t-1];
    for (j=0; j < N_CHANGE_STATS; j++){
      factorval = INPUT_PARAM[j];
      if (hattr == factorval) CHANGE_STAT[j] += s;
      if (tattr == factorval) CHANGE_STAT[j] += s;
    }
  })
}

/*****************
 stat: nonzero
*****************/
WtD_CHANGESTAT_FN(d_nonzero){
  int i;
  
  CHANGE_STAT[0] = 0.0;
  EXEC_THROUGH_TOGGLES(i,{
        CHANGE_STAT[0] += (weights[i]!=0) - (CURWT!=0);
  })
}


/********************  changestats:   S    ***********/

/*****************
 stat: nsumlogfactorial
 Turns a Poisson-reference or geometric-reference ERGM into a Conway-Maxwell-Poisson distribution
*****************/
WtD_CHANGESTAT_FN(d_nsumlogfactorial){
  int i;
  
  CHANGE_STAT[0] = 0.0;
  EXEC_THROUGH_TOGGLES(i,{
      CHANGE_STAT[0] -= lgamma1p(weights[i])-lgamma1p(CURWT);
  })
}

/*****************
 stat: sum
*****************/
WtD_CHANGESTAT_FN(d_sum){
  int i;
  
  CHANGE_STAT[0] = 0.0;
  EXEC_THROUGH_TOGGLES(i,{
      CHANGE_STAT[0] += weights[i]-CURWT;
  })
}

/********************  changestats:   T    ***********/

/*****************
 stat: transitiveweights(_max)
*****************/

WtD_FROM_S_FN(d_transitiveweights_max)

WtS_CHANGESTAT_FN(s_transitiveweights_max){ 
  Edge e1, e2;
  Vertex h, t, node3;
  
  CHANGE_STAT[0]=0;
  for (h=1; h <= N_NODES; h++){
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

WtS_CHANGESTAT_FN(s_transitiveweights_sum){ 
  Edge e1, e2;
  Vertex h, t, node3;
  
  CHANGE_STAT[0]=0;
  for (h=1; h <= N_NODES; h++){
    EXEC_THROUGH_FOUTEDGES(h, e1, t, {
      double path_strength = 0;
      EXEC_THROUGH_INEDGES(t, e2, node3, { 
	path_strength += fmin(GETWT(h,node3),GETWT(node3,t));
	})
      CHANGE_STAT[0] += fmin(path_strength, GETWT(h,t));
      })
  }
}

/*****************
 stat: transitiveweights (threshold)
*****************/

WtD_FROM_S_FN(d_transitiveweights_threshold)

WtS_CHANGESTAT_FN(s_transitiveweights_threshold){ 
  Edge e1, e2;
  Vertex h, t, node3;
  
  CHANGE_STAT[0]=0;
  for (h=1; h <= N_NODES; h++){
    EXEC_THROUGH_FOUTEDGES(h, e1, t, {
	if(GETWT(h,t)<=INPUT_ATTRIB[0]) break;
	EXEC_THROUGH_INEDGES(t, e2, node3, { 
	    if(GETWT(h,node3)>INPUT_ATTRIB[0] && GETWT(node3,t)>INPUT_ATTRIB[0]){
	      CHANGE_STAT[0]++; 
	      break;
	    }
	  })
	  })
      }
}

