#include "wtchangestats.h"

/********************  changestats:   A    ***********/

/*****************                       
 stat: absdiff(_nonzero)
*****************/
WtD_CHANGESTAT_FN(d_absdiff_nonzero){ 
  double p = INPUT_ATTRIB[0];
  
  EXEC_THROUGH_TOGGLES({
      if(p==1.0){
	CHANGE_STAT[0] += fabs(INPUT_ATTRIB[TAIL] - INPUT_ATTRIB[HEAD])*((NEWWT!=0)-(OLDWT!=0));
      } else {
	CHANGE_STAT[0] += pow(fabs(INPUT_ATTRIB[TAIL] - INPUT_ATTRIB[HEAD]), p)*((NEWWT!=0)-(OLDWT!=0));
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
	CHANGE_STAT[0] += fabs(INPUT_ATTRIB[TAIL] - INPUT_ATTRIB[HEAD])*(NEWWT-OLDWT);
      } else {
	CHANGE_STAT[0] += pow(fabs(INPUT_ATTRIB[TAIL] - INPUT_ATTRIB[HEAD]), p)*(NEWWT-OLDWT);
      }
    });
}

/*****************
 stat: absdiffcat(_nonzero)
*****************/
WtD_CHANGESTAT_FN(d_absdiffcat_nonzero){ 
  double change, absdiff, NAsubstitute, tailval, headval;
  Vertex ninputs;
  int j;
  
  ninputs = N_INPUT_PARAMS - N_NODES;
  NAsubstitute = INPUT_PARAM[ninputs-1];
  
  EXEC_THROUGH_TOGGLES({
      change = (NEWWT!=0)-(OLDWT!=0);
      tailval = INPUT_ATTRIB[TAIL-1];
      headval = INPUT_ATTRIB[HEAD-1];
      if (tailval == NAsubstitute ||  headval == NAsubstitute) absdiff = NAsubstitute;
      else absdiff = fabs(tailval - headval);
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
  double change, absdiff, NAsubstitute, tailval, headval;
  Vertex ninputs;
  int j;
  
  ninputs = N_INPUT_PARAMS - N_NODES;
  NAsubstitute = INPUT_PARAM[ninputs-1];
  
  EXEC_THROUGH_TOGGLES({
      change = NEWWT-OLDWT;
      tailval = INPUT_ATTRIB[TAIL-1];
      headval = INPUT_ATTRIB[HEAD-1];
      if (tailval == NAsubstitute ||  headval == NAsubstitute) absdiff = NAsubstitute;
      else absdiff = fabs(tailval - headval);
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
 stat: cyclicalweights(_max)
*****************/

WtD_CHANGESTAT_FN(d_cyclicalweights_max){ 
  CHANGE_STAT[0]=0;

  EXEC_THROUGH_TOGGLES({
      // (TAIL,HEAD) as the focus dyad
      // This means that the strongest 2-path doesn't change.
      double best_path = 0;
      EXEC_THROUGH_OUTEDGES(HEAD, e1, k, yhk, {
	  // If the first leg of the 2-path is less than the best
	  // path so far, we can skip the rest of this and save a
	  // GETWT() call.
	  if(yhk > best_path)
	    best_path = fmax(best_path, fmin(GETWT(k,TAIL), yhk));
	});
      CHANGE_STAT[0] += fmin(best_path, NEWWT) - fmin(best_path, OLDWT);

      // (TAIL,HEAD) as the first link in the 2-path
      // This means that only the strongest 2-path may change.
      EXEC_THROUGH_INEDGES(TAIL, e1, j, yjt, {
	  if(j==HEAD) continue;
	  
	  double old_best_path = 0;
	  double new_best_path = 0;

	  EXEC_THROUGH_INEDGES(j, e2, k, ykj, {
	      double old_ytk = (k==HEAD) ? OLDWT : GETWT(TAIL, k);
	      double new_ytk = (k==HEAD) ? NEWWT : old_ytk; 

	      old_best_path = fmax(old_best_path, fmin(old_ytk, ykj));
	      new_best_path = fmax(new_best_path, fmin(new_ytk, ykj));
	    });
	  CHANGE_STAT[0] += fmin(new_best_path, yjt) - fmin(old_best_path, yjt);
	});

      // (TAIL,HEAD) as the second link of the 2-path
      // This means that only the strongest 2-path may change.
      EXEC_THROUGH_OUTEDGES(HEAD, e1, i, yhi, {
	  if(i==TAIL) continue;
	  
	  double old_best_path = 0;
	  double new_best_path = 0;

	  EXEC_THROUGH_OUTEDGES(i, e2, k, yik, {
	      double old_ykh = (k==TAIL) ? OLDWT : GETWT(k,HEAD);
	      double new_ykh = (k==TAIL) ? NEWWT : old_ykh; 

	      old_best_path = fmax(old_best_path, fmin(old_ykh, yik));
	      new_best_path = fmax(new_best_path, fmin(new_ykh, yik));
	    });
	  CHANGE_STAT[0] += fmin(new_best_path, yhi) - fmin(old_best_path, yhi);
	});
    });
}

WtS_CHANGESTAT_FN(s_cyclicalweights_max){ 
  CHANGE_STAT[0]=0;
  for (Vertex tail=1; tail <= N_NODES; tail++){
    EXEC_THROUGH_FOUTEDGES(tail, e1, head, yth, {
      double best_path = 0;
      EXEC_THROUGH_OUTEDGES(head, e2, node3, yh3, { 
	  // If the second leg of the 2-path is less than the best
	  // path so far, we can skip the rest of this and save a
	  // GETWT() call.
	  if(yh3 > best_path)
	    best_path = fmax(best_path, fmin(GETWT(node3,tail),yh3));
	});
      CHANGE_STAT[0] += fmin(best_path, yth);
      });
  }
}

/*****************
 stat: cyclicalweights(_sum)
*****************/

WtD_FROM_S_FN(d_cyclicalweights_sum)

WtS_CHANGESTAT_FN(s_cyclicalweights_sum){ 
  CHANGE_STAT[0]=0;
  for (Vertex tail=1; tail <= N_NODES; tail++){
    EXEC_THROUGH_FOUTEDGES(tail, e1, head, yth, {
      double path_strength = 0;
      EXEC_THROUGH_OUTEDGES(head, e2, node3, yh3, { 
	path_strength += fmin(GETWT(node3,tail),yh3);
	})
      CHANGE_STAT[0] += fmin(path_strength, yth);
      })
  }
}

/*****************
 stat: cyclicalweights(_threshold)
*****************/

WtD_FROM_S_FN(d_cyclicalweights_threshold)

WtS_CHANGESTAT_FN(s_cyclicalweights_threshold){ 
  CHANGE_STAT[0]=0;
  for (Vertex tail=1; tail <= N_NODES; tail++){
    EXEC_THROUGH_FOUTEDGES(tail, e1, head, yth, {
	if(yth<=INPUT_ATTRIB[0]) break;
	EXEC_THROUGH_OUTEDGES(head, e2, node3, yh3, { 
	    if(yh3>INPUT_ATTRIB[0] && GETWT(node3,tail)>INPUT_ATTRIB[0]){
	      CHANGE_STAT[0]++; 
	      break;
	    }
	  })
	  })
      }
}

/********************  changestats:   G    ***********/

/*****************
 changestat: d_edgecov(_nonzero)
*****************/
WtD_CHANGESTAT_FN(d_edgecov_nonzero) {
  Vertex nrow, noffset;
  
  noffset = BIPARTITE;
  if(noffset > 0){
    nrow = noffset;
  }else{
    nrow = N_NODES;
  }
  
  EXEC_THROUGH_TOGGLES({
      double val = INPUT_ATTRIB[(HEAD-1-noffset)*nrow+(TAIL-1)];  
      CHANGE_STAT[0] += val*((NEWWT!=0)-(OLDWT!=0));
    });
}

/*****************
 changestat: d_edgecov(_sum)
*****************/
WtD_CHANGESTAT_FN(d_edgecov_sum) {
  Vertex nrow, noffset;
  
  noffset = BIPARTITE;
  if(noffset > 0){
    nrow = noffset;
  }else{
    nrow = N_NODES;
  }
  
  EXEC_THROUGH_TOGGLES({
      double val = INPUT_ATTRIB[(HEAD-1-noffset)*nrow+(TAIL-1)];  
      CHANGE_STAT[0] += val*(NEWWT-OLDWT);
    });
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
 stat: mutual (product a.k.a. covariance)
*****************/
WtD_CHANGESTAT_FN(d_mutual_wt_product){
  EXEC_THROUGH_TOGGLES({
      double htweight = GETWT(HEAD,TAIL);
      CHANGE_STAT[0] += (NEWWT*htweight) - (OLDWT*htweight);
    });
}

/*****************
 stat: mutual (geometric mean)
*****************/
WtD_CHANGESTAT_FN(d_mutual_wt_geom_mean){
  EXEC_THROUGH_TOGGLES({
      double htweight = GETWT(HEAD,TAIL);
      CHANGE_STAT[0] += sqrt(NEWWT*htweight) - sqrt(OLDWT*htweight);
    });
}

/*****************
 stat: mutual (minimum)
*****************/
WtD_CHANGESTAT_FN(d_mutual_wt_min){
  EXEC_THROUGH_TOGGLES({
      double htweight = GETWT(HEAD,TAIL);
      CHANGE_STAT[0] += fmin(NEWWT,htweight) - fmin(OLDWT,htweight);
    });
}


/*****************
 stat: mutual (-|yij-yji|)
*****************/
WtD_CHANGESTAT_FN(d_mutual_wt_nabsdiff){
  EXEC_THROUGH_TOGGLES({
      double htweight = GETWT(HEAD,TAIL);
      CHANGE_STAT[0] -= fabs(NEWWT-htweight) - fabs(OLDWT-htweight);
  });
}

/********************  changestats:   N    ***********/

/*****************
 stat: nodecov (nonzero)
*****************/
WtD_CHANGESTAT_FN(d_nodecov_nonzero){ 
  EXEC_THROUGH_TOGGLES({
      CHANGE_STAT[0] += (INPUT_ATTRIB[TAIL-1] + INPUT_ATTRIB[HEAD-1])*((NEWWT!=0)-(OLDWT!=0));
  });
}

/*****************
 stat: node covar 
*****************/
WtD_CHANGESTAT_FN(d_nodecovar){
  EXEC_THROUGH_TOGGLES({
      for(Vertex j=1; j<=N_NODES; j++){
	if(j==HEAD || j==TAIL) continue;
	double ytj = GETWT(TAIL,j);
	CHANGE_STAT[0] += (NEWWT*ytj) - (OLDWT*ytj);
      }

      for(Vertex i=1; i<=N_NODES; i++){
	if(i==TAIL || i==HEAD) continue;
	double yih = GETWT(i,HEAD);
	CHANGE_STAT[0] += (NEWWT*yih) - (OLDWT*yih);
      }
    });
}

/*****************
 stat: nodecov (sum)
*****************/
WtD_CHANGESTAT_FN(d_nodecov_sum){ 
  EXEC_THROUGH_TOGGLES({
      CHANGE_STAT[0] += (INPUT_ATTRIB[TAIL-1] + INPUT_ATTRIB[HEAD-1])*(NEWWT-OLDWT);
  });
}

/*****************
 stat: nodeicov (nonzero)
*****************/
WtD_CHANGESTAT_FN(d_nodeicov_nonzero){ 
  EXEC_THROUGH_TOGGLES({
      CHANGE_STAT[0] += INPUT_ATTRIB[HEAD-1]*((NEWWT!=0)-(OLDWT!=0));
  });
}

/*****************
 stat: nodeicov (sum)
*****************/
WtD_CHANGESTAT_FN(d_nodeicov_sum){ 
  EXEC_THROUGH_TOGGLES({
      CHANGE_STAT[0] += INPUT_ATTRIB[HEAD-1]*(NEWWT-OLDWT);
  });
}

/*****************
 stat: nodeocov (nonzero)
*****************/
WtD_CHANGESTAT_FN(d_nodeocov_nonzero){ 
  EXEC_THROUGH_TOGGLES({
      CHANGE_STAT[0] += INPUT_ATTRIB[TAIL-1]*((NEWWT!=0)-(OLDWT!=0));
  });
}

/*****************
 stat: nodeocov (sum)
*****************/
WtD_CHANGESTAT_FN(d_nodeocov_sum){ 
  EXEC_THROUGH_TOGGLES({
      CHANGE_STAT[0] += INPUT_ATTRIB[TAIL-1]*(NEWWT-OLDWT);
  });
}

/*****************
 stat: nodefactor (nonzero)
*****************/
WtD_CHANGESTAT_FN(d_nodefactor_nonzero){ 
  double s, factorval;
  int j, tailattr, headattr;
  
  
  EXEC_THROUGH_TOGGLES({
      s = (NEWWT!=0) - (OLDWT!=0);
      tailattr = INPUT_ATTRIB[TAIL-1];
      headattr = INPUT_ATTRIB[HEAD-1];
      for (j=0; j < N_CHANGE_STATS; j++){
	factorval = INPUT_PARAM[j];
	if (tailattr == factorval) CHANGE_STAT[j] += s;
	if (headattr == factorval) CHANGE_STAT[j] += s;
      }
    });
}

/*****************
 stat: nodefactor (sum)
*****************/
WtD_CHANGESTAT_FN(d_nodefactor_sum){ 
  double s, factorval;
  int j, tailattr, headattr;
  
  
  EXEC_THROUGH_TOGGLES({
    s = NEWWT - OLDWT;
    tailattr = INPUT_ATTRIB[TAIL-1];
    headattr = INPUT_ATTRIB[HEAD-1];
    for (j=0; j < N_CHANGE_STATS; j++){
      factorval = INPUT_PARAM[j];
      if (tailattr == factorval) CHANGE_STAT[j] += s;
      if (headattr == factorval) CHANGE_STAT[j] += s;
    }
  });
}

/*****************
 stat: node i[n] covar 
*****************/
WtD_CHANGESTAT_FN(d_nodeicovar){
  EXEC_THROUGH_TOGGLES({
      for(Vertex i=1; i<=N_NODES; i++){
	if(i==TAIL || i==HEAD) continue;
	double yih = GETWT(i,HEAD);
	CHANGE_STAT[0] += (NEWWT*yih) - (OLDWT*yih);
      }
    });
}

/*****************
 stat: node i[n] sq[uare]r[oo]t covar[iance] 
*****************/
WtD_CHANGESTAT_FN(d_nodeisqrtcovar){
  EXEC_THROUGH_TOGGLES({
      double sqrtdiff = (sqrt(NEWWT)-sqrt(OLDWT))/(N_NODES-2);
      EXEC_THROUGH_INEDGES(HEAD, e, i, yih, {
	  if(i!=TAIL) 
	    CHANGE_STAT[0] += sqrtdiff*sqrt(yih);
	});
    });
}

/*****************
 stat: nodeifactor (nonzero)
*****************/
WtD_CHANGESTAT_FN(d_nodeifactor_nonzero){ 
  double s, factorval;
  int j, headattr;
  
  EXEC_THROUGH_TOGGLES({
      s = (NEWWT!=0) - (OLDWT!=0);
      headattr = INPUT_ATTRIB[TAIL-1];
      for (j=0; j < N_CHANGE_STATS; j++){
	factorval = INPUT_PARAM[j];
	if (headattr == factorval) CHANGE_STAT[j] += s;
      }
    });
}

/*****************
 stat: nodeifactor (sum)
*****************/
WtD_CHANGESTAT_FN(d_nodeifactor_sum){ 
  double s, factorval;
  int j, headattr;
  
  EXEC_THROUGH_TOGGLES({
    s = NEWWT - OLDWT;
    headattr = INPUT_ATTRIB[TAIL-1];
    for (j=0; j < N_CHANGE_STATS; j++){
      factorval = INPUT_PARAM[j];
      if (headattr == factorval) CHANGE_STAT[j] += s;
    }
  });
}

/*****************
 stat: node o[ut] covar 
*****************/
WtD_CHANGESTAT_FN(d_nodeocovar){
  EXEC_THROUGH_TOGGLES({
      for(Vertex j=1; j<=N_NODES; j++){
	if(j==HEAD || j==TAIL) continue;
	double ytj = GETWT(TAIL,j);
	CHANGE_STAT[0] += (NEWWT*ytj) - (OLDWT*ytj);
      }
    });
}

/*****************
 stat: node o[ut] sq[uare]r[oo]t covar 
*****************/
WtD_CHANGESTAT_FN(d_nodeosqrtcovar){
  EXEC_THROUGH_TOGGLES({
      double sqrtdiff = (sqrt(NEWWT)-sqrt(OLDWT))/(N_NODES-2);
      EXEC_THROUGH_OUTEDGES(TAIL, e, j, ytj, {
	  if(j!=HEAD) 
	    CHANGE_STAT[0] += sqrtdiff*sqrt(ytj);
	});
    });
}

/*****************
 stat: nodeofactor (nonzero)
*****************/
WtD_CHANGESTAT_FN(d_nodeofactor_nonzero){ 
  double s, factorval;
  int j, tailattr;
  
  
  EXEC_THROUGH_TOGGLES({
      s = (NEWWT!=0) - (OLDWT!=0);
      tailattr = INPUT_ATTRIB[TAIL-1];
      for (j=0; j < N_CHANGE_STATS; j++){
	factorval = INPUT_PARAM[j];
	if (tailattr == factorval) CHANGE_STAT[j] += s;
      }
    });
}

/*****************
 stat: nodeofactor (sum)
*****************/
WtD_CHANGESTAT_FN(d_nodeofactor_sum){ 
  double s, factorval;
  int j, tailattr;
  
  
  EXEC_THROUGH_TOGGLES({
    s = NEWWT - OLDWT;
    tailattr = INPUT_ATTRIB[TAIL-1];
    for (j=0; j < N_CHANGE_STATS; j++){
      factorval = INPUT_PARAM[j];
      if (tailattr == factorval) CHANGE_STAT[j] += s;
    }
  });
}

/*****************
 stat: node sq[uare]r[oo]t covar[iance] "centered"
*****************/

WtD_CHANGESTAT_FN(d_nodesqrtcovar_centered){
  // Compute sum(sqrt(y)) (directed) or twice that (undirected). We can update it for each toggle.
  double ssq = 0;
  if(DIRECTED){
    for(Vertex i=1; i<=N_NODES; i++){
      EXEC_THROUGH_EDGES(i, e, j, yij, {
	  ssq += sqrt(yij);
	});
    }
  }else{
    for(Vertex i=1; i<=N_NODES; i++){
      EXEC_THROUGH_FOUTEDGES(i, e, j, yij, {
	  ssq += sqrt(yij);
	});
    }
    ssq*=2;
  }

  EXEC_THROUGH_TOGGLES({
      double change = 0;
      double sqrtdiff = sqrt(NEWWT)-sqrt(OLDWT);
      double new_ssq = ssq + sqrtdiff*(DIRECTED? 1 : 2);
      
      EXEC_THROUGH_EDGES(TAIL, e, j, ytj, {
	  if(j!=HEAD) change += sqrtdiff*sqrt(ytj);
	});
      
      EXEC_THROUGH_EDGES(HEAD, e, i, yih, {
	  if(i!=TAIL) change += sqrtdiff*sqrt(yih);
	});

      CHANGE_STAT[0] += change/(N_NODES-2);
      CHANGE_STAT[0] -= (new_ssq*new_ssq-ssq*ssq) / (N_NODES*(N_NODES-1)) / 2;

      ssq = new_ssq;
    });
}

WtS_CHANGESTAT_FN(s_nodesqrtcovar_centered){
  CHANGE_STAT[0] = 0;
  double ssq = 0;
  for(Vertex i=1; i<=N_NODES; i++){
    EXEC_THROUGH_EDGES(i, e1, j, yij, {
	double sqrtyij = sqrt(yij);
	ssq += sqrtyij;
	EXEC_THROUGH_EDGES(i, e2, k, yik, {
	    if(k>=j) break; // sqrt(yij)*sqrt(yik)==sqrt(yik)*sqrt(yij)
	    CHANGE_STAT[0] += sqrtyij*sqrt(yik);
	  });
      });
  }
  CHANGE_STAT[0] /= N_NODES-2;
  CHANGE_STAT[0] -= ssq*ssq / (N_NODES*(N_NODES-1)) / 2; // (The /2 is because of k<j.)
}

/*****************
 stat: node sq[uare]r[oo]t covar[iance]
*****************/
WtD_CHANGESTAT_FN(d_nodesqrtcovar){
  EXEC_THROUGH_TOGGLES({
      double change = 0;
      double sqrtdiff = sqrt(NEWWT)-sqrt(OLDWT);
      
      EXEC_THROUGH_EDGES(TAIL, e, j, ytj, {
	  if(j!=HEAD) change += sqrtdiff*sqrt(ytj);
	});
      
      EXEC_THROUGH_EDGES(HEAD, e, i, yih, {
	  if(i!=TAIL) change += sqrtdiff*sqrt(yih);
	});

      CHANGE_STAT[0] += change/(N_NODES-2);
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

WtD_CHANGESTAT_FN(d_transitiveweights_max){ 
  CHANGE_STAT[0]=0;

  EXEC_THROUGH_TOGGLES({
      // (TAIL,HEAD) as the focus dyad
      // This means that the strongest 2-path doesn't change.
      double best_path = 0;
      EXEC_THROUGH_INEDGES(HEAD, e1, k, ykh, {
	  // If the second leg of the 2-path is less than the best
	  // path so far, we can skip the rest of this and save a
	  // GETWT() call.
	  if(ykh > best_path)
	    best_path = fmax(best_path, fmin(GETWT(TAIL,k), ykh));
	});
      CHANGE_STAT[0] += fmin(best_path, NEWWT) - fmin(best_path, OLDWT);

      // (TAIL,HEAD) as the first link of the 2-path
      // This means that only the strongest 2-path may change.
      EXEC_THROUGH_OUTEDGES(TAIL, e1, j, ytj, {
	  if(j==HEAD) continue;
	  
	  double old_best_path = 0;
	  double new_best_path = 0;

	  EXEC_THROUGH_INEDGES(j, e2, k, ykj, {
	      double old_ytk = (k==HEAD) ? OLDWT : GETWT(TAIL, k);
	      double new_ytk = (k==HEAD) ? NEWWT : old_ytk; 

	      old_best_path = fmax(old_best_path, fmin(old_ytk, ykj));
	      new_best_path = fmax(new_best_path, fmin(new_ytk, ykj));
	    });
	  CHANGE_STAT[0] += fmin(new_best_path, ytj) - fmin(old_best_path, ytj);
	});

      // (TAIL,HEAD) as the second link of the 2-path
      // This means that only the strongest 2-path may change.
      EXEC_THROUGH_INEDGES(HEAD, e1, i, yih, {
	  if(i==TAIL) continue;
	  
	  double old_best_path = 0;
	  double new_best_path = 0;

	  EXEC_THROUGH_OUTEDGES(i, e2, k, yik, {
	      double old_ykh = (k==TAIL) ? OLDWT : GETWT(k,HEAD);
	      double new_ykh = (k==TAIL) ? NEWWT : old_ykh; 

	      old_best_path = fmax(old_best_path, fmin(old_ykh, yik));
	      new_best_path = fmax(new_best_path, fmin(new_ykh, yik));
	    });
	  CHANGE_STAT[0] += fmin(new_best_path, yih) - fmin(old_best_path, yih);
	});
    });
}

WtS_CHANGESTAT_FN(s_transitiveweights_max){ 
  CHANGE_STAT[0]=0;
  for (Vertex tail=1; tail <= N_NODES; tail++){
    EXEC_THROUGH_FOUTEDGES(tail, e1, head, yth, {
      double best_path = 0;
      EXEC_THROUGH_INEDGES(head, e2, node3, y3h, {
	  // If the second leg of the 2-path is less than the best
	  // path so far, we can skip the rest of this and save a
	  // GETWT() call.
	  if(y3h > best_path)
	    best_path = fmax(best_path, fmin(GETWT(tail,node3), y3h));
	});
      CHANGE_STAT[0] += fmin(best_path, yth);
      });
  }
}

/*****************
 stat: transitiveweights(_sum)
*****************/

WtD_FROM_S_FN(d_transitiveweights_sum)

WtS_CHANGESTAT_FN(s_transitiveweights_sum){ 
  CHANGE_STAT[0]=0;
  for (Vertex tail=1; tail <= N_NODES; tail++){
    EXEC_THROUGH_FOUTEDGES(tail, e1, head, yth, {
      double path_strength = 0;
      EXEC_THROUGH_INEDGES(head, e2, node3, y3h, { 
	path_strength += fmin(GETWT(tail,node3),y3h);
	})
      CHANGE_STAT[0] += fmin(path_strength, yth);
      })
  }
}

/*****************
 stat: transitiveweights (threshold)
*****************/

WtD_FROM_S_FN(d_transitiveweights_threshold)

WtS_CHANGESTAT_FN(s_transitiveweights_threshold){ 
  CHANGE_STAT[0]=0;
  for (Vertex tail=1; tail <= N_NODES; tail++){
    EXEC_THROUGH_FOUTEDGES(tail, e1, head, yth, {
	if(yth<=INPUT_ATTRIB[0]) continue;
	EXEC_THROUGH_INEDGES(head, e2, node3, y3h, { 
	    if(y3h>INPUT_ATTRIB[0] && GETWT(tail,node3)>INPUT_ATTRIB[0]){
	      CHANGE_STAT[0]++; 
	      break;
	    }
	  });
      });
  }
}

