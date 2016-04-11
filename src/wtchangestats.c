/*  File src/wtchangestats.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2013 Statnet Commons
 */
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

/*****************
 stat: atmost
*****************/
WtD_CHANGESTAT_FN(d_atmost){
  EXEC_THROUGH_TOGGLES({
      CHANGE_STAT[0] += (NEWWT<=INPUT_ATTRIB[0]) - (OLDWT<=INPUT_ATTRIB[0]);
    });
}

/********************  changestats:   B    ***********/

/*****************
 stat: b2cov (nonzero)
*****************/
WtD_CHANGESTAT_FN(d_b2cov_nonzero){ 
  EXEC_THROUGH_TOGGLES({
      CHANGE_STAT[0] += INPUT_ATTRIB[HEAD-BIPARTITE-1]*((NEWWT!=0)-(OLDWT!=0));
  });
}

/*****************
 stat: b2cov (sum)
*****************/
WtD_CHANGESTAT_FN(d_b2cov_sum){ 
  EXEC_THROUGH_TOGGLES({
      CHANGE_STAT[0] += INPUT_ATTRIB[HEAD-BIPARTITE-1]*(NEWWT-OLDWT);
  });
}

/*****************
 stat: b2factor (nonzero)
*****************/
WtD_CHANGESTAT_FN(d_b2factor_nonzero){ 
  double s, factorval;
  int j, headattr;
  
  EXEC_THROUGH_TOGGLES({
      s = (NEWWT!=0) - (OLDWT!=0);
      headattr = INPUT_ATTRIB[HEAD-BIPARTITE-1];
      for (j=0; j < N_CHANGE_STATS; j++){
	factorval = INPUT_PARAM[j];
	if (headattr == factorval) CHANGE_STAT[j] += s;
      }
    });
}

/*****************
 stat: b2factor (sum)
*****************/
WtD_CHANGESTAT_FN(d_b2factor_sum){ 
  double s, factorval;
  int j, headattr;
  
  EXEC_THROUGH_TOGGLES({
    s = NEWWT - OLDWT;
    headattr = INPUT_ATTRIB[HEAD-BIPARTITE-1];
    for (j=0; j < N_CHANGE_STATS; j++){
      factorval = INPUT_PARAM[j];
      if (headattr == factorval) CHANGE_STAT[j] += s;
    }
  });
}


/********************  changestats:   C    ***********/

/*****************
 stat: cyclicalweights
*****************/

WtD_CHANGESTAT_FN(d_cyclicalweights){ 
  unsigned int path = INPUT_ATTRIB[0], combine = INPUT_ATTRIB[1], compare =  INPUT_ATTRIB[2];
  CHANGE_STAT[0]=0;

  EXEC_THROUGH_TOGGLES({
      // (TAIL,HEAD) as the focus dyad
      // This means that the strongest 2-path doesn't change.
      double two_paths = 0;
      EXEC_THROUGH_OUTEDGES(HEAD, e1, k, yhk, {
	  double two_path;

	  switch(path){
	  case 1: two_path = fmin(GETWT(k,TAIL), yhk); break; // min
	  case 2: two_path = sqrt(GETWT(k,TAIL) * yhk); break; // geomean
	  default: // never reached, but prevents a warning
	    two_path = 0;
	  }
	  
	  switch(combine){
	  case 1: two_paths = fmax(two_paths, two_path); break; // max
	  case 2: two_paths += two_path; break; // sum
	  }
	});

      switch(compare){
      case 1: CHANGE_STAT[0] += fmin(two_paths, NEWWT) - fmin(two_paths, OLDWT); break; // min
      case 2: CHANGE_STAT[0] += sqrt(two_paths * NEWWT) - sqrt(two_paths * OLDWT); break; // geomean
      }

      // (TAIL,HEAD) as the first link in the 2-path
      // This means that only the strongest 2-path may change.
      EXEC_THROUGH_INEDGES(TAIL, e1, j, yjt, {
	  if(j==HEAD) continue;
	  
	  double old_two_paths = 0;
	  double new_two_paths = 0;

	  EXEC_THROUGH_INEDGES(j, e2, k, ykj, {
	      double old_ytk = (k==HEAD) ? OLDWT : GETWT(TAIL, k);
	      double new_ytk = (k==HEAD) ? NEWWT : old_ytk; 
	      double old_two_path;
	      double new_two_path;

	      switch(path){
	      case 1: // min
		old_two_path = fmin(old_ytk, ykj);
		new_two_path = fmin(new_ytk, ykj);
		break;
	      case 2: // geomean
		old_two_path = sqrt(old_ytk * ykj);
		new_two_path = sqrt(new_ytk * ykj);
		break;
	      default: // never reached, but prevents a warning
		old_two_path = 0;
		new_two_path = 0;
	      }
	  
	      switch(combine){
	      case 1:
		old_two_paths = fmax(old_two_paths, old_two_path);// max
		new_two_paths = fmax(new_two_paths, new_two_path);
		break; // max
	      case 2:
		old_two_paths += old_two_path;
		new_two_paths += new_two_path;
		break; // sum
	      }

	    });

	  switch(compare){
	  case 1: CHANGE_STAT[0] += fmin(new_two_paths, yjt) - fmin(old_two_paths, yjt); break; // min
	  case 2: CHANGE_STAT[0] += sqrt(new_two_paths * yjt) - sqrt(old_two_paths * yjt); break; // geomean
	  }
	});

      // (TAIL,HEAD) as the second link of the 2-path
      // This means that only the strongest 2-path may change.
      EXEC_THROUGH_OUTEDGES(HEAD, e1, i, yhi, {
	  if(i==TAIL) continue;
	  
	  double old_two_paths = 0;
	  double new_two_paths = 0;

	  EXEC_THROUGH_OUTEDGES(i, e2, k, yik, {
	      double old_ykh = (k==TAIL) ? OLDWT : GETWT(k,HEAD);
	      double new_ykh = (k==TAIL) ? NEWWT : old_ykh; 
	      double old_two_path;
	      double new_two_path;

	      switch(path){
	      case 1: // min
		old_two_path = fmin(old_ykh, yik);
		new_two_path = fmin(new_ykh, yik);
		break;
	      case 2: // geomean
		old_two_path = sqrt(old_ykh * yik);
		new_two_path = sqrt(new_ykh * yik);
		break;
	      default: // never reached, but prevents a warning
		old_two_path = 0;
		new_two_path = 0;
	      }
	  
	      switch(combine){
	      case 1:
		old_two_paths = fmax(old_two_paths, old_two_path);// max
		new_two_paths = fmax(new_two_paths, new_two_path);
		break; // max
	      case 2:
		old_two_paths += old_two_path;
		new_two_paths += new_two_path;
		break; // sum
	      }
	    });
	  
	  switch(compare){
	  case 1: CHANGE_STAT[0] += fmin(new_two_paths, yhi) - fmin(old_two_paths, yhi); break; // min
	  case 2: CHANGE_STAT[0] += sqrt(new_two_paths * yhi) - sqrt(old_two_paths * yhi); break; // geomean
	  }
	});
    });
}

WtS_CHANGESTAT_FN(s_cyclicalweights){ 
  unsigned int path = INPUT_ATTRIB[0], combine = INPUT_ATTRIB[1], compare = INPUT_ATTRIB[2];
  /* unsigned int threshold = INPUT_ATTRIB[3]; */

  CHANGE_STAT[0]=0;
  for (Vertex tail=1; tail <= N_NODES; tail++){
    EXEC_THROUGH_FOUTEDGES(tail, e1, head, yth, {
      double two_paths = 0;
      EXEC_THROUGH_OUTEDGES(head, e2, node3, yh3, { 
	  double two_path;
	  
	  switch(path){
	  case 1: two_path = fmin(GETWT(node3,tail), yh3); break; // min
	  case 2: two_path = sqrt(GETWT(node3,tail) * yh3); break; // geomean
	  default: // never reached, but prevents a warning
	    two_path = 0;
	  }
	  
	  switch(combine){
	  case 1: two_paths = fmax(two_paths, two_path); break; // max
	  case 2: two_paths += two_path; break; // sum
	  }

	});
    
      switch(compare){
      case 1: CHANGE_STAT[0] += fmin(two_paths, yth); break; // min
      case 2: CHANGE_STAT[0] += sqrt(two_paths * yth); break; // geomean
      }
      
      });
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

/********************  changestats:   L    ***********/

/********************  changestats:   M    ***********/

/*****************
 changestat: d_mix_nonzero
 This appears to be the version of nodemix used for 
 bipartite networks (only)
*****************/
WtD_CHANGESTAT_FN(d_mix_nonzero){
  int matchvaltail, matchvalhead;
  int j, nstats;

  nstats = N_CHANGE_STATS;

  /* *** don't forget tail -> head */    
  EXEC_THROUGH_TOGGLES({
      double s = (NEWWT!=0) - (OLDWT!=0);
      matchvaltail = INPUT_PARAM[TAIL-1+2*nstats];
      matchvalhead = INPUT_PARAM[HEAD-1+2*nstats];
      for (j=0; j<nstats; j++) {
	if(matchvaltail==INPUT_PARAM[j] && matchvalhead==INPUT_PARAM[nstats+j]) {
	  CHANGE_STAT[j] += s;
	}
      }
    });
}

/*****************
 changestat: d_mix_sum
 This appears to be the version of nodemix used for 
 bipartite networks (only)
*****************/
WtD_CHANGESTAT_FN(d_mix_sum){
  int matchvaltail, matchvalhead;
  int j, nstats;

  nstats = N_CHANGE_STATS;

  /* *** don't forget tail -> head */    
  EXEC_THROUGH_TOGGLES({
      double s = NEWWT - OLDWT;
      matchvaltail = INPUT_PARAM[TAIL-1+2*nstats];
      matchvalhead = INPUT_PARAM[HEAD-1+2*nstats];
      for (j=0; j<nstats; j++) {
	if(matchvaltail==INPUT_PARAM[j] && matchvalhead==INPUT_PARAM[nstats+j]) {
	  CHANGE_STAT[j] += s;
	}
      }
    });
}

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
      headattr = INPUT_ATTRIB[HEAD-1];
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
    headattr = INPUT_ATTRIB[HEAD-1];
    for (j=0; j < N_CHANGE_STATS; j++){
      factorval = INPUT_PARAM[j];
      if (headattr == factorval) CHANGE_STAT[j] += s;
    }
  });
}

/*****************
 changestat: d_nodematch_nonzero
*****************/
WtD_CHANGESTAT_FN(d_nodematch_nonzero) { 
  Vertex ninputs = N_INPUT_PARAMS - N_NODES;

  /* *** don't forget tail -> head */    
  EXEC_THROUGH_TOGGLES({
      double matchval = INPUT_PARAM[TAIL+ninputs-1];
      if (matchval == INPUT_PARAM[HEAD+ninputs-1]) { /* We have a match! */
	double s = (NEWWT!=0) - (OLDWT!=0);
	if (ninputs==0) {/* diff=F in network statistic specification */
	  CHANGE_STAT[0] += s;
	} else { /* diff=T */
	  for (unsigned int j=0; j<ninputs; j++) {
	    if (matchval == INPUT_PARAM[j]) 
	      CHANGE_STAT[j] += s;
	  }
	}
      }
    });
}

/*****************
 changestat: d_nodematch_sum
*****************/
WtD_CHANGESTAT_FN(d_nodematch_sum) { 
  Vertex ninputs = N_INPUT_PARAMS - N_NODES;

  /* *** don't forget tail -> head */    
  EXEC_THROUGH_TOGGLES({
      double matchval = INPUT_PARAM[TAIL+ninputs-1];
      if (matchval == INPUT_PARAM[HEAD+ninputs-1]) { /* We have a match! */
	double s = NEWWT - OLDWT;
	if (ninputs==0) {/* diff=F in network statistic specification */
	  CHANGE_STAT[0] += s;
	} else { /* diff=T */
	  for (unsigned int j=0; j<ninputs; j++) {
	    if (matchval == INPUT_PARAM[j]) 
	      CHANGE_STAT[j] += s;
	  }
	}
      }
    });
}

/*****************
 changestat: d_nodemix_nonzero
 Update mixing matrix, non-bipartite networks only 
 (but see also d_mix_nonzero)
*****************/
WtD_CHANGESTAT_FN(d_nodemix_nonzero) {
  int ninputs = N_INPUT_PARAMS - N_NODES, ninputs2 = ninputs/2;
  double rtype, ctype, tmp;

  /* *** don't forget tail -> head */    
  EXEC_THROUGH_TOGGLES({
      double change = (NEWWT!=0) - (OLDWT!=0);
      /*Find the node covariate values (types) for the tail and head*/
      rtype=INPUT_PARAM[TAIL+ninputs-1];
      ctype=INPUT_PARAM[HEAD+ninputs-1];
      if (!DIRECTED && rtype > ctype)  {
        tmp = rtype; rtype = ctype; ctype = tmp; /* swap rtype, ctype */
      }
      /*Find the right statistic to update */
      for(unsigned int j=0; j<ninputs2; j++){
        if((INPUT_PARAM[j] == rtype) && (INPUT_PARAM[j+ninputs2] == ctype)){
          CHANGE_STAT[j] += change;
          j = ninputs2; /* leave the for loop */
        }
      } 
    });
}

/*****************
 changestat: d_nodemix_sum
 Update mixing matrix, non-bipartite networks only 
 (but see also d_mix_sum)
*****************/
WtD_CHANGESTAT_FN(d_nodemix_sum) {
  int ninputs = N_INPUT_PARAMS - N_NODES, ninputs2 = ninputs/2;
  double rtype, ctype, tmp;

  /* *** don't forget tail -> head */    
  EXEC_THROUGH_TOGGLES({
      double change = NEWWT - OLDWT;
      /*Find the node covariate values (types) for the tail and head*/
      rtype=INPUT_PARAM[TAIL+ninputs-1];
      ctype=INPUT_PARAM[HEAD+ninputs-1];
      if (!DIRECTED && rtype > ctype)  {
        tmp = rtype; rtype = ctype; ctype = tmp; /* swap rtype, ctype */
      }
      /*Find the right statistic to update */
      for(unsigned int j=0; j<ninputs2; j++){
        if((INPUT_PARAM[j] == rtype) && (INPUT_PARAM[j+ninputs2] == ctype)){
          CHANGE_STAT[j] += change;
          j = ninputs2; /* leave the for loop */
        }
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
	  ssq += sqrt(yij); j=j; /* j=j is silly, just to prevent compiler warnings */
	});
    }
  }else{
    for(Vertex i=1; i<=N_NODES; i++){
      EXEC_THROUGH_FOUTEDGES(i, e, j, yij, {
	  ssq += sqrt(yij); j=j; /* j=j is silly, just to prevent compiler warnings */
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

/********************  changestats:   S    ***********/

/*****************
 stat: smallerthan
*****************/
WtD_CHANGESTAT_FN(d_smallerthan){
  EXEC_THROUGH_TOGGLES({
      CHANGE_STAT[0] += (NEWWT<INPUT_ATTRIB[0]) - (OLDWT<INPUT_ATTRIB[0]);
  });
}

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
 stat: transitiveweights
*****************/

WtD_CHANGESTAT_FN(d_transitiveweights){ 
  CHANGE_STAT[0]=0;
  unsigned int path = INPUT_ATTRIB[0], combine = INPUT_ATTRIB[1], compare =  INPUT_ATTRIB[2];
  
  EXEC_THROUGH_TOGGLES({
      // (TAIL,HEAD) as the focus dyad
      // This means that the combined 2-path value doesn't change.
      double two_paths = 0;
      EXEC_THROUGH_INEDGES(HEAD, e1, k, ykh, {
	  double two_path;

	  switch(path){
	  case 1: two_path = fmin(GETWT(TAIL,k), ykh); break; // min
	  case 2: two_path = sqrt(GETWT(TAIL,k) * ykh); break; // geomean
	  default: // never reached, but prevents a warning
	    two_path = 0;    
	  }
	  
	  switch(combine){
	  case 1: two_paths = fmax(two_paths, two_path); break; // max
	  case 2: two_paths += two_path; break; // sum
	  }

	});

      switch(compare){
      case 1: CHANGE_STAT[0] += fmin(two_paths, NEWWT) - fmin(two_paths, OLDWT); break; // min
      case 2: CHANGE_STAT[0] += sqrt(two_paths * NEWWT) - sqrt(two_paths * OLDWT); break; // geomean
      }

      // (TAIL,HEAD) as the first link of the 2-path
      // This means that only the combined 2-path value may change.
      EXEC_THROUGH_OUTEDGES(TAIL, e1, j, ytj, {
	  if(j==HEAD) continue;
	  
	  double old_two_paths = 0;
	  double new_two_paths = 0;

	  EXEC_THROUGH_INEDGES(j, e2, k, ykj, {
	      double old_ytk = (k==HEAD) ? OLDWT : GETWT(TAIL, k);
	      double new_ytk = (k==HEAD) ? NEWWT : old_ytk; 
	      double old_two_path;
	      double new_two_path;

	      switch(path){
	      case 1: // min
		old_two_path = fmin(old_ytk, ykj);
		new_two_path = fmin(new_ytk, ykj);
		break;
	      case 2: // geomean
		old_two_path = sqrt(old_ytk * ykj);
		new_two_path = sqrt(new_ytk * ykj);
		break;
	      default: // never reached, but prevents a warning
		old_two_path = 0;
		new_two_path = 0;
	      }
	  
	      switch(combine){
	      case 1:
		old_two_paths = fmax(old_two_paths, old_two_path);// max
		new_two_paths = fmax(new_two_paths, new_two_path);
		break; // max
	      case 2:
		old_two_paths += old_two_path;
		new_two_paths += new_two_path;
		break; // sum
	      }
	      
	    });

	  switch(compare){
	  case 1: CHANGE_STAT[0] += fmin(new_two_paths, ytj) - fmin(old_two_paths, ytj); break; // min
	  case 2: CHANGE_STAT[0] += sqrt(new_two_paths * ytj) - sqrt(old_two_paths * ytj); break; // geomean
	  }
	});

      // (TAIL,HEAD) as the second link of the 2-path
      // This means that only the combined 2-path value may change.
      EXEC_THROUGH_INEDGES(HEAD, e1, i, yih, {
	  if(i==TAIL) continue;
	  
	  double old_two_paths = 0;
	  double new_two_paths = 0;

	  EXEC_THROUGH_OUTEDGES(i, e2, k, yik, {
	      double old_ykh = (k==TAIL) ? OLDWT : GETWT(k,HEAD);
	      double new_ykh = (k==TAIL) ? NEWWT : old_ykh; 
	      double old_two_path;
	      double new_two_path;

	      switch(path){
	      case 1: // min
		old_two_path = fmin(old_ykh, yik);
		new_two_path = fmin(new_ykh, yik);
		break;
	      case 2: // geomean
		old_two_path = sqrt(old_ykh * yik);
		new_two_path = sqrt(new_ykh * yik);
		break;
	      default: // never reached, but prevents a warning
		old_two_path = 0;
		new_two_path = 0;
	      }
	  
	      switch(combine){
	      case 1:
		old_two_paths = fmax(old_two_paths, old_two_path);// max
		new_two_paths = fmax(new_two_paths, new_two_path);
		break; // max
	      case 2:
		old_two_paths += old_two_path;
		new_two_paths += new_two_path;
		break; // sum
	      }
	    });

	  switch(compare){
	  case 1: CHANGE_STAT[0] += fmin(new_two_paths, yih) - fmin(old_two_paths, yih); break; // min
	  case 2: CHANGE_STAT[0] += sqrt(new_two_paths * yih) - sqrt(old_two_paths * yih); break; // geomean
	  }
	});
    });
}

WtS_CHANGESTAT_FN(s_transitiveweights){ 
  unsigned int path = INPUT_ATTRIB[0], combine = INPUT_ATTRIB[1], compare =  INPUT_ATTRIB[2];

  CHANGE_STAT[0]=0;
  for (Vertex tail=1; tail <= N_NODES; tail++){
    EXEC_THROUGH_FOUTEDGES(tail, e1, head, yth, {
      double two_paths = 0;
      EXEC_THROUGH_INEDGES(head, e2, node3, y3h, {
	  double two_path;

	  switch(path){
	  case 1: two_path = fmin(GETWT(tail,node3), y3h); break; // min
	  case 2: two_path = sqrt(GETWT(tail,node3) * y3h); break; // geomean
	  default: // never reached, but prevents a warning
	    two_path = 0;    
	  }
	  
	  switch(combine){
	  case 1: two_paths = fmax(two_paths, two_path); break; // max
	  case 2: two_paths += two_path; break; // sum
	  }

	});

      switch(compare){
      case 1: CHANGE_STAT[0] += fmin(two_paths, yth); break; // min
      case 2: CHANGE_STAT[0] += sqrt(two_paths * yth); break; // geomean
      }

      });
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

