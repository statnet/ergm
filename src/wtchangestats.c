#include "wtchangestats.h"

/********************  changestats:   A    ***********/

/*****************                       
 stat: absdiff(_nonzero)
*****************/
WtD_CHANGESTAT_FN(d_absdiff_nonzero){ 
  double p = INPUT_ATTRIB[0];
  
  EXEC_THROUGH_TOGGLES({
      if(p==1.0){
	CHANGE_STAT[0] += fabs(INPUT_ATTRIB[HEAD] - INPUT_ATTRIB[TAIL])*((NEWWT!=0)-(OLDWT!=0));
      } else {
	CHANGE_STAT[0] += pow(fabs(INPUT_ATTRIB[HEAD] - INPUT_ATTRIB[TAIL]), p)*((NEWWT!=0)-(OLDWT!=0));
      }
    });
}

/*****************                       
 stat: absdiff(_sum)
*****************/
WtD_CHANGESTAT_FN(d_absdiff_sum){ 
  double p = INPUT_ATTRIB[0];
  
  EXEC_THROUGH_TOGGLES({
      if(p==1.0){
	CHANGE_STAT[0] += fabs(INPUT_ATTRIB[HEAD] - INPUT_ATTRIB[TAIL])*(NEWWT-OLDWT);
      } else {
	CHANGE_STAT[0] += pow(fabs(INPUT_ATTRIB[HEAD] - INPUT_ATTRIB[TAIL]), p)*(NEWWT-OLDWT);
      }
    });
}

/*****************
 stat: absdiffcat(_nonzero)
*****************/
WtD_CHANGESTAT_FN(d_absdiffcat_nonzero){ 
  double change, absdiff, NAsubstitute, hval, tval;
  Vertex ninputs;
  int j;
  
  ninputs = N_INPUT_PARAMS - N_NODES;
  NAsubstitute = INPUT_PARAM[ninputs-1];
  
  EXEC_THROUGH_TOGGLES({
      change = (NEWWT!=0)-(OLDWT!=0);
      hval = INPUT_ATTRIB[HEAD-1];
      tval = INPUT_ATTRIB[TAIL-1];
      if (hval == NAsubstitute ||  tval == NAsubstitute) absdiff = NAsubstitute;
      else absdiff = fabs(hval - tval);
      if (absdiff>0){
	for (j=0; j<N_CHANGE_STATS; j++){
	  CHANGE_STAT[j] += (absdiff==INPUT_PARAM[j]) ? change : 0.0;
	}
      }
    });
}

/*****************
 stat: absdiffcat(_sum)
*****************/
WtD_CHANGESTAT_FN(d_absdiffcat_sum){ 
  double change, absdiff, NAsubstitute, hval, tval;
  Vertex ninputs;
  int j;
  
  ninputs = N_INPUT_PARAMS - N_NODES;
  NAsubstitute = INPUT_PARAM[ninputs-1];
  
  EXEC_THROUGH_TOGGLES({
      change = NEWWT-OLDWT;
      hval = INPUT_ATTRIB[HEAD-1];
      tval = INPUT_ATTRIB[TAIL-1];
      if (hval == NAsubstitute ||  tval == NAsubstitute) absdiff = NAsubstitute;
      else absdiff = fabs(hval - tval);
      if (absdiff>0){
	for (j=0; j<N_CHANGE_STATS; j++){
	  CHANGE_STAT[j] += (absdiff==INPUT_PARAM[j]) ? change : 0.0;
	}
      }
    });
}


/*****************
 stat: atleast
*****************/
WtD_CHANGESTAT_FN(d_atleast){
  EXEC_THROUGH_TOGGLES({
      CHANGE_STAT[0] += (NEWWT>=INPUT_ATTRIB[0]) - (OLDWT>=INPUT_ATTRIB[0]);
    });
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
  EXEC_THROUGH_TOGGLES({
      CHANGE_STAT[0] += (NEWWT>INPUT_ATTRIB[0]) - (OLDWT>INPUT_ATTRIB[0]);
  });
}

/********************  changestats:   I    ***********/

/*****************
 stat: ininterval
*****************/
WtD_CHANGESTAT_FN(d_ininterval){
  EXEC_THROUGH_TOGGLES({
      CHANGE_STAT[0] += ((INPUT_ATTRIB[2] ? NEWWT>INPUT_ATTRIB[0] : NEWWT>=INPUT_ATTRIB[0]) && (INPUT_ATTRIB[3] ? NEWWT<INPUT_ATTRIB[1] : NEWWT<=INPUT_ATTRIB[1])) - ((INPUT_ATTRIB[2] ? OLDWT>INPUT_ATTRIB[0] : OLDWT>=INPUT_ATTRIB[0]) && (INPUT_ATTRIB[3] ? OLDWT<INPUT_ATTRIB[1] : OLDWT<=INPUT_ATTRIB[1]));
    });
}


/********************  changestats:   M    ***********/

/*****************
 stat: mutual (product a.k.a. correlation)
*****************/
WtD_CHANGESTAT_FN(d_mutual_wt_product){
  EXEC_THROUGH_TOGGLES({
      double thweight = GETWT(TAIL,HEAD);
      CHANGE_STAT[0] += (NEWWT*thweight) - (OLDWT*thweight);
    });
}

/*****************
 stat: mutual (geometric mean)
*****************/
WtD_CHANGESTAT_FN(d_mutual_wt_geom_mean){
  EXEC_THROUGH_TOGGLES({
      double thweight = GETWT(TAIL,HEAD);
      CHANGE_STAT[0] += sqrt(NEWWT*thweight) - sqrt(OLDWT*thweight);
    });
}

/*****************
 stat: mutual (minimum)
*****************/
WtD_CHANGESTAT_FN(d_mutual_wt_min){
  EXEC_THROUGH_TOGGLES({
      double thweight = GETWT(TAIL,HEAD);
      CHANGE_STAT[0] += fmin(NEWWT,thweight) - fmin(OLDWT,thweight);
    });
}


/*****************
 stat: mutual (-|yij-yji|)
*****************/
WtD_CHANGESTAT_FN(d_mutual_wt_nabsdiff){
  EXEC_THROUGH_TOGGLES({
      double thweight = GETWT(TAIL,HEAD);
      CHANGE_STAT[0] -= fabs(NEWWT-thweight) - fabs(OLDWT-thweight);
  });
}

/********************  changestats:   N    ***********/

/*****************
 stat: nodecov (nonzero)
*****************/
WtD_CHANGESTAT_FN(d_nodecov_nonzero){ 
  EXEC_THROUGH_TOGGLES({
      CHANGE_STAT[0] += (INPUT_ATTRIB[HEAD-1] + INPUT_ATTRIB[TAIL-1])*((NEWWT!=0)-(OLDWT!=0));
  });
}

/*****************
 stat: d_nodecov
*****************/
WtD_CHANGESTAT_FN(d_nodecov_sum){ 
  EXEC_THROUGH_TOGGLES({
      CHANGE_STAT[0] += (INPUT_ATTRIB[HEAD-1] + INPUT_ATTRIB[TAIL-1])*(NEWWT-OLDWT);
  });
}

/*****************
 stat: nodefactor (nonzero)
*****************/
WtD_CHANGESTAT_FN(d_nodefactor_nonzero){ 
  double s, factorval;
  int j, hattr, tattr;
  
  
  EXEC_THROUGH_TOGGLES({
      s = (NEWWT!=0) - (OLDWT!=0);
      hattr = INPUT_ATTRIB[HEAD-1];
      tattr = INPUT_ATTRIB[TAIL-1];
      for (j=0; j < N_CHANGE_STATS; j++){
	factorval = INPUT_PARAM[j];
	if (hattr == factorval) CHANGE_STAT[j] += s;
	if (tattr == factorval) CHANGE_STAT[j] += s;
      }
    });
}

/*****************
 stat: nodefactor (sum)
*****************/
WtD_CHANGESTAT_FN(d_nodefactor_sum){ 
  double s, factorval;
  int j, hattr, tattr;
  
  
  EXEC_THROUGH_TOGGLES({
    s = NEWWT - OLDWT;
    hattr = INPUT_ATTRIB[HEAD-1];
    tattr = INPUT_ATTRIB[TAIL-1];
    for (j=0; j < N_CHANGE_STATS; j++){
      factorval = INPUT_PARAM[j];
      if (hattr == factorval) CHANGE_STAT[j] += s;
      if (tattr == factorval) CHANGE_STAT[j] += s;
    }
  });
}

/*****************
 stat: node i[n] corr 
*****************/
WtD_CHANGESTAT_FN(d_nodeicorr){
  EXEC_THROUGH_TOGGLES({
      for(Vertex i=1; i<=N_NODES; i++){
	if(i==HEAD) continue;
	double yit = GETWT(i,TAIL);
	CHANGE_STAT[0] += (NEWWT*yit) - (OLDWT*yit);
      }
    });
}

/*****************
 stat: nodeifactor (nonzero)
*****************/
WtD_CHANGESTAT_FN(d_nodeifactor_nonzero){ 
  double s, factorval;
  int j, tattr;
  
  
  EXEC_THROUGH_TOGGLES({
      s = (NEWWT!=0) - (OLDWT!=0);
      tattr = INPUT_ATTRIB[HEAD-1];
      for (j=0; j < N_CHANGE_STATS; j++){
	factorval = INPUT_PARAM[j];
	if (tattr == factorval) CHANGE_STAT[j] += s;
      }
    });
}

/*****************
 stat: nodeifactor (sum)
*****************/
WtD_CHANGESTAT_FN(d_nodeifactor_sum){ 
  double s, factorval;
  int j, tattr;
  
  
  EXEC_THROUGH_TOGGLES({
    s = NEWWT - OLDWT;
    tattr = INPUT_ATTRIB[HEAD-1];
    for (j=0; j < N_CHANGE_STATS; j++){
      factorval = INPUT_PARAM[j];
      if (tattr == factorval) CHANGE_STAT[j] += s;
    }
  });
}

/*****************
 stat: node o[ut] corr 
*****************/
WtD_CHANGESTAT_FN(d_nodeocorr){
  EXEC_THROUGH_TOGGLES({
      for(Vertex j=1; j<=N_NODES; j++){
	if(j==TAIL) continue;
	double yhj = GETWT(HEAD,j);
	CHANGE_STAT[0] += (NEWWT*yhj) - (OLDWT*yhj);
      }
    });
}

/*****************
 stat: nodeofactor (nonzero)
*****************/
WtD_CHANGESTAT_FN(d_nodeofactor_nonzero){ 
  double s, factorval;
  int j, hattr;
  
  
  EXEC_THROUGH_TOGGLES({
      s = (NEWWT!=0) - (OLDWT!=0);
      hattr = INPUT_ATTRIB[HEAD-1];
      for (j=0; j < N_CHANGE_STATS; j++){
	factorval = INPUT_PARAM[j];
	if (hattr == factorval) CHANGE_STAT[j] += s;
      }
    });
}

/*****************
 stat: nodeofactor (sum)
*****************/
WtD_CHANGESTAT_FN(d_nodeofactor_sum){ 
  double s, factorval;
  int j, hattr;
  
  
  EXEC_THROUGH_TOGGLES({
    s = NEWWT - OLDWT;
    hattr = INPUT_ATTRIB[HEAD-1];
    for (j=0; j < N_CHANGE_STATS; j++){
      factorval = INPUT_PARAM[j];
      if (hattr == factorval) CHANGE_STAT[j] += s;
    }
  });
}


/*****************
 stat: nonzero
*****************/
WtD_CHANGESTAT_FN(d_nonzero){
  EXEC_THROUGH_TOGGLES({
        CHANGE_STAT[0] += (NEWWT!=0) - (OLDWT!=0);
  });
}

/*****************
 stat: nsumlogfactorial
 Turns a Poisson-reference or geometric-reference ERGM into a Conway-Maxwell-Poisson distribution
*****************/
WtD_CHANGESTAT_FN(d_nsumlogfactorial){
  EXEC_THROUGH_TOGGLES({
      CHANGE_STAT[0] -= lgamma1p(NEWWT)-lgamma1p(OLDWT);
  });
}

/********************  changestats:   S    ***********/

/*****************
 stat: sum
*****************/
WtD_CHANGESTAT_FN(d_sum){
  EXEC_THROUGH_TOGGLES({
      CHANGE_STAT[0] += NEWWT-OLDWT;
  });
}

/*****************
 stat: sum(with power)
*****************/
WtD_CHANGESTAT_FN(d_sum_pow){
  double p = INPUT_ATTRIB[0];
  
  EXEC_THROUGH_TOGGLES({
      CHANGE_STAT[0] += pow(NEWWT,p)-pow(OLDWT,p);
    });
}

/********************  changestats:   T    ***********/

/*****************
 stat: transitiveweights(_max)
*****************/

/* WtD_CHANGESTAT_FN(d_transitiveweights_max){ */
/*   Edge e1, e2; */
/*   Vertex node3; */
/*   EXEC_THROUGH_TOGGLES({ */
/*       /\* Changed dyad as the focus dyad. *\/ */
/*       double best_path = 0; */
/*       EXEC_THROUGH_INEDGES(TAIL, e1, node3, {  */
/* 	best_path = fmax(best_path, fmin(GETWT(HEAD,node3),GETWT(node3,TAIL))); */
/* 	}) */
/*       CHANGE_STAT[0] += fmin(best_path, NEWWT) - fmin(best_path, OLDWT); */

/*       /\* Changed dyad as a part of a two-path.  */
/* 	 A dyad (i,j) is potentially affected by (h,t) iff: */
/* 	 y(i,j)>0 & (h=i & (t,j)>0 | t=j & (i,h)>0). */
/*       *\/ */
/*       /\* For all ties (h=i,j)>0, *\/ */
/*       EXEC_THROUGH_OUTEDGES(HEAD, e1, node3, { */
/* 	  double ijwt = GETWT(TAIL,node3); */
/* 	  /\* If (t,j)>0), (h,t) can affect (i,j). *\/ */
/* 	  if(ijwt>0){ */
/* 	    best_path = 0; */
/* 	    EXEC_THROUGH_INEDGES(TAIL, e1, node3, {  */
/* 		best_path = fmax(best_path, fmin(GETWT(HEAD,node3),GETWT(node3,TAIL))); */
/* 	      }) */
/* 	      CHANGE_STAT[0] += fmin(best_path, NEWWT) - fmin(best_path, OLDWT); */
/* 	  } */
/* 	} */
/*   }); */
/* } */

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

