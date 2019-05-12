/*  File src/wtchangestats.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2019 Statnet Commons
 */
#include "wtchangestats.h"

/********************  changestats:   A    ***********/

/*****************                       
 stat: absdiff(_nonzero)
*****************/
WtC_CHANGESTAT_FN(c_absdiff_nonzero){ 
  double p = INPUT_ATTRIB[0];
  
  ZERO_ALL_CHANGESTATS();
      if(p==1.0){
	CHANGE_STAT[0] += fabs(INPUT_ATTRIB[tail] - INPUT_ATTRIB[head])*((weight!=0)-(GETWT(tail,head)!=0));
      } else {
	CHANGE_STAT[0] += pow(fabs(INPUT_ATTRIB[tail] - INPUT_ATTRIB[head]), p)*((weight!=0)-(GETWT(tail,head)!=0));
      }
}

/*****************                       
 stat: absdiff(_sum)
*****************/
WtC_CHANGESTAT_FN(c_absdiff_sum){ 
  double p = INPUT_ATTRIB[0];
  
  ZERO_ALL_CHANGESTATS();
      if(p==1.0){
	CHANGE_STAT[0] += fabs(INPUT_ATTRIB[tail] - INPUT_ATTRIB[head])*(weight-GETWT(tail,head));
      } else {
	CHANGE_STAT[0] += pow(fabs(INPUT_ATTRIB[tail] - INPUT_ATTRIB[head]), p)*(weight-GETWT(tail,head));
      }
}

/*****************
 stat: absdiffcat(_nonzero)
*****************/
WtC_CHANGESTAT_FN(c_absdiffcat_nonzero){ 
  double change, absdiff, NAsubstitute, tailval, headval;
  Vertex ninputs;
  int j;
  
  ninputs = N_INPUT_PARAMS - N_NODES;
  NAsubstitute = INPUT_PARAM[ninputs-1];
  
  ZERO_ALL_CHANGESTATS();
      change = (weight!=0)-(GETWT(tail,head)!=0);
      tailval = INPUT_ATTRIB[tail-1];
      headval = INPUT_ATTRIB[head-1];
      if (tailval == NAsubstitute ||  headval == NAsubstitute) absdiff = NAsubstitute;
      else absdiff = fabs(tailval - headval);
      if (absdiff>0){
	for (j=0; j<N_CHANGE_STATS; j++){
	  CHANGE_STAT[j] += (absdiff==INPUT_PARAM[j]) ? change : 0.0;
	}
      }
}

/*****************
 stat: absdiffcat(_sum)
*****************/
WtC_CHANGESTAT_FN(c_absdiffcat_sum){ 
  double change, absdiff, NAsubstitute, tailval, headval;
  Vertex ninputs;
  int j;
  
  ninputs = N_INPUT_PARAMS - N_NODES;
  NAsubstitute = INPUT_PARAM[ninputs-1];
  
  ZERO_ALL_CHANGESTATS();
      change = weight-GETWT(tail,head);
      tailval = INPUT_ATTRIB[tail-1];
      headval = INPUT_ATTRIB[head-1];
      if (tailval == NAsubstitute ||  headval == NAsubstitute) absdiff = NAsubstitute;
      else absdiff = fabs(tailval - headval);
      if (absdiff>0){
	for (j=0; j<N_CHANGE_STATS; j++){
	  CHANGE_STAT[j] += (absdiff==INPUT_PARAM[j]) ? change : 0.0;
	}
      }
}


/*****************
 stat: atleast
*****************/
WtC_CHANGESTAT_FN(c_atleast){
  ZERO_ALL_CHANGESTATS();
  double oldwt = GETWT(tail,head);
  for(unsigned int i=0; i<N_CHANGE_STATS; i++)
    CHANGE_STAT[i] += (weight>=INPUT_ATTRIB[i]) - (oldwt>=INPUT_ATTRIB[i]);
}

/*****************
 stat: atmost
*****************/
WtC_CHANGESTAT_FN(c_atmost){
  ZERO_ALL_CHANGESTATS();
  double oldwt = GETWT(tail,head);
  for(unsigned int i=0; i<N_CHANGE_STATS; i++)
    CHANGE_STAT[i] += (weight<=INPUT_ATTRIB[i]) - (oldwt<=INPUT_ATTRIB[i]);
}

/********************  changestats:   B    ***********/

/*****************
 stat: b2cov (nonzero)
*****************/
WtC_CHANGESTAT_FN(c_b2cov_nonzero){ 
  unsigned int oshift = N_INPUT_PARAMS / N_CHANGE_STATS;
  ZERO_ALL_CHANGESTATS();
  double oldwt = GETWT(tail,head);
  for(unsigned int j=0, o=0; j<N_CHANGE_STATS; j++, o+=oshift)
    CHANGE_STAT[j] += INPUT_ATTRIB[head-BIPARTITE+o-1]*((weight!=0)-(oldwt!=0));
}

/*****************
 stat: b2cov (sum)
*****************/
WtC_CHANGESTAT_FN(c_b2cov_sum){ 
  unsigned int oshift = N_INPUT_PARAMS / N_CHANGE_STATS;
  ZERO_ALL_CHANGESTATS();
  double oldwt = GETWT(tail,head);
  for(unsigned int j=0, o=0; j<N_CHANGE_STATS; j++, o+=oshift)
    CHANGE_STAT[j] += INPUT_ATTRIB[head-BIPARTITE+o-1]*(weight-oldwt);
}

/*****************
 stat: b2factor (nonzero)
*****************/
WtC_CHANGESTAT_FN(c_b2factor_nonzero){
  ZERO_ALL_CHANGESTATS();
  double oldwt = GETWT(tail,head);
  double s = (weight!=0) - (oldwt!=0);
  int headpos = INPUT_ATTRIB[head-1-BIPARTITE];
  if(headpos != -1) CHANGE_STAT[headpos] += s;
}

/*****************
 stat: b2factor (sum)
*****************/
WtC_CHANGESTAT_FN(c_b2factor_sum){
  ZERO_ALL_CHANGESTATS();
  double oldwt = GETWT(tail,head);
  double s = weight - oldwt;
  int headpos = INPUT_ATTRIB[head-1-BIPARTITE];
  if(headpos != -1) CHANGE_STAT[headpos] += s;
}

/********************  changestats:   C    ***********/

/*****************
 stat: cyclicalweights
*****************/

WtC_CHANGESTAT_FN(c_cyclicalweights){ 
  unsigned int path = INPUT_ATTRIB[0], combine = INPUT_ATTRIB[1], compare =  INPUT_ATTRIB[2];
  CHANGE_STAT[0]=0;

  ZERO_ALL_CHANGESTATS();
      // (tail,head) as the focus dyad
      // This means that the strongest 2-path doesn't change.
      double two_paths = 0;
      EXEC_THROUGH_OUTEDGES(head, e1, k, yhk, {
	  double two_path;

	  switch(path){
	  case 1: two_path = fmin(GETWT(k,tail), yhk); break; // min
	  case 2: two_path = sqrt(GETWT(k,tail) * yhk); break; // geomean
	  default: // never reached, but prevents a warning
	    two_path = 0;
	  }
	  
	  switch(combine){
	  case 1: two_paths = fmax(two_paths, two_path); break; // max
	  case 2: two_paths += two_path; break; // sum
	  }
	});

      switch(compare){
      case 1: CHANGE_STAT[0] += fmin(two_paths, weight) - fmin(two_paths, GETWT(tail,head)); break; // min
      case 2: CHANGE_STAT[0] += sqrt(two_paths * weight) - sqrt(two_paths * GETWT(tail,head)); break; // geomean
      }

      // (tail,head) as the first link in the 2-path
      // This means that only the strongest 2-path may change.
      EXEC_THROUGH_INEDGES(tail, e1, j, yjt, {
	  if(j==head) continue;
	  
	  double old_two_paths = 0;
	  double new_two_paths = 0;

	  EXEC_THROUGH_INEDGES(j, e2, k, ykj, {
	      double old_ytk = (k==head) ? GETWT(tail,head) : GETWT(tail, k);
	      double new_ytk = (k==head) ? weight : old_ytk; 
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

      // (tail,head) as the second link of the 2-path
      // This means that only the strongest 2-path may change.
      EXEC_THROUGH_OUTEDGES(head, e1, i, yhi, {
	  if(i==tail) continue;
	  
	  double old_two_paths = 0;
	  double new_two_paths = 0;

	  EXEC_THROUGH_OUTEDGES(i, e2, k, yik, {
	      double old_ykh = (k==tail) ? GETWT(tail,head) : GETWT(k,head);
	      double new_ykh = (k==tail) ? weight : old_ykh; 
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
WtC_CHANGESTAT_FN(c_edgecov_nonzero) {
  Vertex nrow, noffset;
  
  noffset = BIPARTITE;
  if(noffset > 0){
    nrow = noffset;
  }else{
    nrow = N_NODES;
  }
  
  ZERO_ALL_CHANGESTATS();
      double val = INPUT_ATTRIB[(head-1-noffset)*nrow+(tail-1)];  
      CHANGE_STAT[0] += val*((weight!=0)-(GETWT(tail,head)!=0));
}
/*****************
 changestat: d_edgecov(_sum)
*****************/
WtC_CHANGESTAT_FN(c_edgecov_sum) {
  Vertex nrow, noffset;
  
  noffset = BIPARTITE;
  if(noffset > 0){
    nrow = noffset;
  }else{
    nrow = N_NODES;
  }
  
  ZERO_ALL_CHANGESTATS();
      double val = INPUT_ATTRIB[(head-1-noffset)*nrow+(tail-1)];  
      CHANGE_STAT[0] += val*(weight-GETWT(tail,head));
}
/********************  changestats:   D    ***********/

/*****************                       
 changestat: d_diff(_nonzero)
*****************/
WtC_CHANGESTAT_FN(c_diff_nonzero) { 
  double p = INPUT_PARAM[0], *x = INPUT_PARAM+2;
  int mul = INPUT_PARAM[1], sign_code = INPUT_PARAM[2];

  /* *** don't forget tail -> head */
  ZERO_ALL_CHANGESTATS();
    double change = (x[tail] - x[head])*mul;
    switch(sign_code){
    case 1: // identity
      break;
    case 2: // abs
      change = fabs(change);
      break;
    case 3: // positive only
      change = change<0 ? 0 : change;
      break;
    case 4: // negative only
      change = change>0 ? 0 : change;
      break;
    default:
      error("Invalid sign action code passed to d_diff_nonzero.");
      break;
    }

    if(p==0.0){ // Special case: take the sign of the difference instead.
      change = sign(change);
    }else if(p!=1.0){
      change = pow(change, p);
    }
    
    CHANGE_STAT[0] += change*((weight!=0)-(GETWT(tail,head)!=0));
}

/*****************                       
 changestat: d_diff(_sum)
*****************/
WtC_CHANGESTAT_FN(c_diff_sum) { 
  double p = INPUT_PARAM[0], *x = INPUT_PARAM+2;
  int mul = INPUT_PARAM[1], sign_code = INPUT_PARAM[2];

  /* *** don't forget tail -> head */
  ZERO_ALL_CHANGESTATS();
    double change = (x[tail] - x[head])*mul;
    switch(sign_code){
    case 1: // identity
      break;
    case 2: // abs
      change = fabs(change);
      break;
    case 3: // positive only
      change = change<0 ? 0 : change;
      break;
    case 4: // negative only
      change = change>0 ? 0 : change;
      break;
    default:
      error("Invalid sign action code passed to d_diff_sum.");
      break;
    }

    if(p==0.0){ // Special case: take the sign of the difference instead.
      change = sign(change);
    }else if(p!=1.0){
      change = pow(change, p);
    }
    
    CHANGE_STAT[0] += change*(weight-GETWT(tail,head));
}

    
/********************  changestats:   G    ***********/

/*****************
 stat: greaterthan
*****************/
WtC_CHANGESTAT_FN(c_greaterthan){
  ZERO_ALL_CHANGESTATS();
  double oldwt = GETWT(tail,head);
  for(unsigned int i=0; i<N_CHANGE_STATS; i++)
    CHANGE_STAT[i] += (weight>INPUT_ATTRIB[i]) - (oldwt>INPUT_ATTRIB[i]);
}

/********************  changestats:   I    ***********/

/*****************
 stat: ininterval
*****************/
WtC_CHANGESTAT_FN(c_ininterval){
  ZERO_ALL_CHANGESTATS();
      CHANGE_STAT[0] += ((INPUT_ATTRIB[2] ? weight>INPUT_ATTRIB[0] : weight>=INPUT_ATTRIB[0]) && (INPUT_ATTRIB[3] ? weight<INPUT_ATTRIB[1] : weight<=INPUT_ATTRIB[1])) - ((INPUT_ATTRIB[2] ? GETWT(tail,head)>INPUT_ATTRIB[0] : GETWT(tail,head)>=INPUT_ATTRIB[0]) && (INPUT_ATTRIB[3] ? GETWT(tail,head)<INPUT_ATTRIB[1] : GETWT(tail,head)<=INPUT_ATTRIB[1]));
}
/********************  changestats:   L    ***********/

/********************  changestats:   M    ***********/

/*****************
 changestat: d_mix_nonzero
 This appears to be the version of nodemix used for 
 bipartite networks (only)
*****************/
WtC_CHANGESTAT_FN(c_mix_nonzero){
  int matchvaltail, matchvalhead;
  int j, nstats;

  nstats = N_CHANGE_STATS;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS();
      double s = (weight!=0) - (GETWT(tail,head)!=0);
      matchvaltail = INPUT_PARAM[tail-1+2*nstats];
      matchvalhead = INPUT_PARAM[head-1+2*nstats];
      for (j=0; j<nstats; j++) {
	if(matchvaltail==INPUT_PARAM[j] && matchvalhead==INPUT_PARAM[nstats+j]) {
	  CHANGE_STAT[j] += s;
	}
      }
}
/*****************
 changestat: d_mix_sum
 This appears to be the version of nodemix used for 
 bipartite networks (only)
*****************/
WtC_CHANGESTAT_FN(c_mix_sum){
  int matchvaltail, matchvalhead;
  int j, nstats;

  nstats = N_CHANGE_STATS;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS();
      double s = weight - GETWT(tail,head);
      matchvaltail = INPUT_PARAM[tail-1+2*nstats];
      matchvalhead = INPUT_PARAM[head-1+2*nstats];
      for (j=0; j<nstats; j++) {
	if(matchvaltail==INPUT_PARAM[j] && matchvalhead==INPUT_PARAM[nstats+j]) {
	  CHANGE_STAT[j] += s;
	}
      }
}
/*****************
 stat: mutual (product a.k.a. covariance)
*****************/
WtC_CHANGESTAT_FN(c_mutual_wt_product){
  ZERO_ALL_CHANGESTATS();
      double htweight = GETWT(head,tail);
      CHANGE_STAT[0] += (weight*htweight) - (GETWT(tail,head)*htweight);
}
/*****************
 stat: mutual (geometric mean)
*****************/
WtC_CHANGESTAT_FN(c_mutual_wt_geom_mean){
  ZERO_ALL_CHANGESTATS();
      double htweight = GETWT(head,tail);
      CHANGE_STAT[0] += sqrt(weight*htweight) - sqrt(GETWT(tail,head)*htweight);
}
/*****************
 stat: mutual (minimum)
*****************/
WtC_CHANGESTAT_FN(c_mutual_wt_min){
  ZERO_ALL_CHANGESTATS();
      double htweight = GETWT(head,tail);
      CHANGE_STAT[0] += fmin(weight,htweight) - fmin(GETWT(tail,head),htweight);
}

/*****************
 stat: mutual (-|yij-yji|)
*****************/
WtC_CHANGESTAT_FN(c_mutual_wt_nabsdiff){
  ZERO_ALL_CHANGESTATS();
      double htweight = GETWT(head,tail);
      CHANGE_STAT[0] -= fabs(weight-htweight) - fabs(GETWT(tail,head)-htweight);
}

/********************  changestats:   N    ***********/

/*****************
 stat: nodecov (nonzero)
*****************/
WtC_CHANGESTAT_FN(c_nodecov_nonzero){ 
  unsigned int oshift = N_INPUT_PARAMS / N_CHANGE_STATS;
  ZERO_ALL_CHANGESTATS();
  double oldwt = GETWT(tail,head);
  for(unsigned int j=0, o=0; j<N_CHANGE_STATS; j++, o+=oshift)
    CHANGE_STAT[j] += (INPUT_ATTRIB[tail+o-1] + INPUT_ATTRIB[head+o-1])*((weight!=0)-(oldwt!=0));
}

/*****************
 stat: node covar[iance] 
*****************/
WtC_CHANGESTAT_FN(c_nodecovar){
  unsigned int transcode = INPUT_ATTRIB[0], center = INPUT_ATTRIB[1];
  double sum = center?*(double *)STORAGE : 0;
  
  ZERO_ALL_CHANGESTATS();
      double diff = TRANSFORM_DYADVAL(weight,transcode)-TRANSFORM_DYADVAL(GETWT(tail,head),transcode);
      double new_sum = center? sum + diff : 0;
      diff /= (N_NODES-2);
      EXEC_THROUGH_EDGES(tail, e, i, yti, {
	  if(i!=head)
	    CHANGE_STAT[0] += diff*TRANSFORM_DYADVAL(yti,transcode);
	});
      EXEC_THROUGH_EDGES(head, e, i, yih, {
	  if(i!=tail)
	    CHANGE_STAT[0] += diff*TRANSFORM_DYADVAL(yih,transcode);
	});
      if(center){
	CHANGE_STAT[0] += (sum*sum-new_sum*new_sum)/N_DYADS;
	sum = new_sum;
      }
}
WtI_CHANGESTAT_FN(i_nodecovar){
  unsigned int transcode = INPUT_ATTRIB[0], center = INPUT_ATTRIB[1];
  if(center){
    ALLOC_STORAGE(1, double, sum);
    *sum = 0;
    EXEC_THROUGH_NET_EDGES(tail, e1, head, y, {
	*sum+=TRANSFORM_DYADVAL(y, transcode);
	head=head; e1=e1; // Prevent a compiler warning.
      });
  }
}

WtU_CHANGESTAT_FN(u_nodecovar){
  unsigned int transcode = INPUT_ATTRIB[0], center = INPUT_ATTRIB[1];
  if(center){
    GET_STORAGE(double, sum);
    if(tail) *sum += TRANSFORM_DYADVAL(weight,transcode)-TRANSFORM_DYADVAL(GETWT(tail,head),transcode);
  }
}

/*****************
 stat: nodecov (sum)
*****************/
WtC_CHANGESTAT_FN(c_nodecov_sum){ 
  unsigned int oshift = N_INPUT_PARAMS / N_CHANGE_STATS;
  ZERO_ALL_CHANGESTATS();
  double oldwt = GETWT(tail,head);
  for(unsigned int j=0, o=0; j<N_CHANGE_STATS; j++, o+=oshift)
    CHANGE_STAT[j] += (INPUT_ATTRIB[tail+o-1] + INPUT_ATTRIB[head+o-1])*(weight-oldwt);
}

/*****************
 stat: nodeicov (nonzero)
*****************/
WtC_CHANGESTAT_FN(c_nodeicov_nonzero){ 
  unsigned int oshift = N_INPUT_PARAMS / N_CHANGE_STATS;
  ZERO_ALL_CHANGESTATS();
  double oldwt = GETWT(tail,head);
  for(unsigned int j=0, o=0; j<N_CHANGE_STATS; j++, o+=oshift)
    CHANGE_STAT[j] += INPUT_ATTRIB[head+o-1]*((weight!=0)-(oldwt!=0));
}

/*****************
 stat: nodeicov (sum)
*****************/
WtC_CHANGESTAT_FN(c_nodeicov_sum){ 
  unsigned int oshift = N_INPUT_PARAMS / N_CHANGE_STATS;
  ZERO_ALL_CHANGESTATS();
  double oldwt = GETWT(tail,head);
  for(unsigned int j=0, o=0; j<N_CHANGE_STATS; j++, o+=oshift)
    CHANGE_STAT[j] += INPUT_ATTRIB[head+o-1]*(weight-oldwt);
}

/*****************
 stat: nodeocov (nonzero)
*****************/
WtC_CHANGESTAT_FN(c_nodeocov_nonzero){ 
  unsigned int oshift = N_INPUT_PARAMS / N_CHANGE_STATS;
  ZERO_ALL_CHANGESTATS();
  double oldwt = GETWT(tail,head);
  for(unsigned int j=0, o=0; j<N_CHANGE_STATS; j++, o+=oshift)
    CHANGE_STAT[j] += INPUT_ATTRIB[tail+o-1]*((weight!=0)-(oldwt!=0));
}

/*****************
 stat: nodeocov (sum)
*****************/
WtC_CHANGESTAT_FN(c_nodeocov_sum){ 
  unsigned int oshift = N_INPUT_PARAMS / N_CHANGE_STATS;
  ZERO_ALL_CHANGESTATS();
  double oldwt = GETWT(tail,head);
  for(unsigned int j=0, o=0; j<N_CHANGE_STATS; j++, o+=oshift)
    CHANGE_STAT[j] += INPUT_ATTRIB[tail+o-1]*(weight-oldwt);
}

/*****************
 stat: nodefactor (nonzero)
*****************/
WtC_CHANGESTAT_FN(c_nodefactor_nonzero){ 
  double s = (weight!=0) - (GETWT(tail,head)!=0);
  int tailpos = INPUT_ATTRIB[tail-1];
  int headpos = INPUT_ATTRIB[head-1];
  if (tailpos != -1) CHANGE_STAT[tailpos] += s;
  if (headpos != -1) CHANGE_STAT[headpos] += s;
}
/*****************
 stat: nodefactor (sum)
*****************/
WtC_CHANGESTAT_FN(c_nodefactor_sum){ 
  double s = weight - GETWT(tail,head);
  int tailpos = INPUT_ATTRIB[tail-1];
  int headpos = INPUT_ATTRIB[head-1];
  if (tailpos != -1) CHANGE_STAT[tailpos] += s;
  if (headpos != -1) CHANGE_STAT[headpos] += s;
}

/*****************
 changestat: receiver (nonzero)
*****************/
WtC_CHANGESTAT_FN(c_receiver_nonzero){ 
  double s = (weight!=0) - (GETWT(tail,head)!=0);
  unsigned int j=0;
  Vertex deg = (Vertex)INPUT_PARAM[j];
  while((deg != head) && (j < (N_CHANGE_STATS-1))){
    j++;
    deg = (Vertex)INPUT_PARAM[j];
  }
  if(deg==head) CHANGE_STAT[j] += s;
}
/*****************
 changestat: receiver (sum)
*****************/
WtC_CHANGESTAT_FN(c_receiver_sum){ 
  double s = weight - GETWT(tail,head);
  unsigned int j=0;
  Vertex deg = (Vertex)INPUT_PARAM[j];
  while((deg != head) && (j < (N_CHANGE_STATS-1))){
    j++;
    deg = (Vertex)INPUT_PARAM[j];
  }
  if(deg==head) CHANGE_STAT[j] += s;
}

/*****************
 changestat: sender (nonzero)
*****************/
WtC_CHANGESTAT_FN(c_sender_nonzero){ 
  double s = (weight!=0) - (GETWT(tail,head)!=0);
  unsigned int j=0;
  Vertex deg = (Vertex)INPUT_PARAM[j];
  while((deg != tail) && (j < (N_CHANGE_STATS-1))){
    j++;
    deg = (Vertex)INPUT_PARAM[j];
  }
  if(deg==tail) CHANGE_STAT[j] += s;
}
/*****************
 changestat: sender (sum)
*****************/
WtC_CHANGESTAT_FN(c_sender_sum){ 
  double s = weight - GETWT(tail,head);
  unsigned int j=0;
  Vertex deg = (Vertex)INPUT_PARAM[j];
  while((deg != tail) && (j < (N_CHANGE_STATS-1))){
    j++;
    deg = (Vertex)INPUT_PARAM[j];
  }
  if(deg==tail) CHANGE_STAT[j] += s;
}

/*****************
 changestat: sociality (nonzero)
*****************/
WtC_CHANGESTAT_FN(c_sociality_nonzero) { 
  int ninputs, nstats;
  
  ninputs = (int)N_INPUT_PARAMS;
  nstats  = (int)N_CHANGE_STATS;

  if(ninputs>nstats+1){
    /* match on attributes */
    double s = (weight!=0) - (GETWT(tail,head)!=0);
    double tailattr = INPUT_ATTRIB[tail-1+nstats+1]; // +1 for the "guard" value between vertex IDs and attribute vector
    if(tailattr == INPUT_ATTRIB[head-1+nstats+1]){
      unsigned int j=0;
      Vertex deg = (Vertex)INPUT_PARAM[j];
      while(deg != tail && j < nstats){
	j++;
	deg = (Vertex)INPUT_PARAM[j];
      }
      if(j < nstats) CHANGE_STAT[j] += s;
      j=0;
      deg = (Vertex)INPUT_PARAM[j];
      while(deg != head && j < nstats){
	j++;
	deg = (Vertex)INPUT_PARAM[j];
      }
      if(j < nstats) CHANGE_STAT[j] += s;
    }
  }else{
    /* *** don't forget tail -> head */    
    double s = (weight!=0) - (GETWT(tail,head)!=0);
    unsigned int j=0;
    Vertex deg = (Vertex)INPUT_PARAM[j];
    while(deg != tail && j < nstats){
      j++;
      deg = (Vertex)INPUT_PARAM[j];
    }
    if(j < nstats) CHANGE_STAT[j] += s;
    j=0;
    deg = (Vertex)INPUT_PARAM[j];
    while(deg != head && j < nstats){
      j++;
      deg = (Vertex)INPUT_PARAM[j];
    }
    if(j < nstats) CHANGE_STAT[j] += s;
  }
}
/*****************
 changestat: sociality (sum)
*****************/
WtC_CHANGESTAT_FN(c_sociality_sum) { 
  int ninputs, nstats;
  
  ninputs = (int)N_INPUT_PARAMS;
  nstats  = (int)N_CHANGE_STATS;

  if(ninputs>nstats+1){
    /* match on attributes */
    double s = weight - GETWT(tail,head);
    double tailattr = INPUT_ATTRIB[tail-1+nstats+1]; // +1 for the "guard" value between vertex IDs and attribute vector
    if(tailattr == INPUT_ATTRIB[head-1+nstats+1]){
      unsigned int j=0;
      Vertex deg = (Vertex)INPUT_PARAM[j];
      while(deg != tail && j < nstats){
	j++;
	deg = (Vertex)INPUT_PARAM[j];
      }
      if(j < nstats) CHANGE_STAT[j] += s;
      j=0;
      deg = (Vertex)INPUT_PARAM[j];
      while(deg != head && j < nstats){
	j++;
	deg = (Vertex)INPUT_PARAM[j];
      }
      if(j < nstats) CHANGE_STAT[j] += s;
    }
  }else{
    /* *** don't forget tail -> head */    
    double s = weight - GETWT(tail,head);
    unsigned int j=0;
    Vertex deg = (Vertex)INPUT_PARAM[j];
    while(deg != tail && j < nstats){
      j++;
      deg = (Vertex)INPUT_PARAM[j];
    }
    if(j < nstats) CHANGE_STAT[j] += s;
    j=0;
    deg = (Vertex)INPUT_PARAM[j];
    while(deg != head && j < nstats){
      j++;
      deg = (Vertex)INPUT_PARAM[j];
    }
    if(j < nstats) CHANGE_STAT[j] += s;
  }
}

/*****************
 stat: node i[n] covar 
*****************/
WtC_CHANGESTAT_FN(c_nodeicovar){
  unsigned int transcode = INPUT_ATTRIB[0], center = INPUT_ATTRIB[1];
  double sum = center?*(double *)STORAGE : 0;
  
  ZERO_ALL_CHANGESTATS();
      double diff = TRANSFORM_DYADVAL(weight,transcode)-TRANSFORM_DYADVAL(GETWT(tail,head),transcode);
      double new_sum = center? sum + diff : 0;
      diff /= (N_NODES-2);
      EXEC_THROUGH_INEDGES(head, e, i, yih, {
	  if(i!=tail)
	    CHANGE_STAT[0] += 2*diff*TRANSFORM_DYADVAL(yih,transcode);
	});
      if(center){
	CHANGE_STAT[0] += (sum*sum-new_sum*new_sum)/N_DYADS;
	sum = new_sum;
      }
}
WtI_CHANGESTAT_FN(i_nodeicovar){
  unsigned int transcode = INPUT_ATTRIB[0], center = INPUT_ATTRIB[1];
  if(center){
    ALLOC_STORAGE(1, double, sum);
    *sum = 0;
    EXEC_THROUGH_NET_EDGES(tail, e1, head, y, {
	*sum+=TRANSFORM_DYADVAL(y, transcode);
	head=head; e1=e1; // Prevent a compiler warning.
      });
  }
}

WtU_CHANGESTAT_FN(u_nodeicovar){
  unsigned int transcode = INPUT_ATTRIB[0], center = INPUT_ATTRIB[1];
  if(center){
    GET_STORAGE(double, sum);
    
    if(tail) *sum += TRANSFORM_DYADVAL(weight,transcode)-TRANSFORM_DYADVAL(GETWT(tail,head),transcode);
  }
}

/*****************
 stat: nodeifactor (nonzero)
*****************/
WtC_CHANGESTAT_FN(c_nodeifactor_nonzero){ 
  double s = (weight!=0) - (GETWT(tail,head)!=0);
  int headpos = INPUT_ATTRIB[head-1];
  if (headpos != -1) CHANGE_STAT[headpos] += s;
}
/*****************
 stat: nodeifactor (sum)
*****************/
WtC_CHANGESTAT_FN(c_nodeifactor_sum){ 
  double s = weight - GETWT(tail,head);
  int headpos = INPUT_ATTRIB[head-1];
  if (headpos != -1) CHANGE_STAT[headpos] += s;
}

/*****************
 changestat: d_nodematch_nonzero
*****************/
WtC_CHANGESTAT_FN(c_nodematch_nonzero) { 
  Vertex ninputs = N_INPUT_PARAMS - N_NODES;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS();
      double matchval = INPUT_PARAM[tail+ninputs-1];
      if (matchval == INPUT_PARAM[head+ninputs-1]) { /* We have a match! */
	double s = (weight!=0) - (GETWT(tail,head)!=0);
	if (ninputs==0) {/* diff=F in network statistic specification */
	  CHANGE_STAT[0] += s;
	} else { /* diff=T */
	  for (unsigned int j=0; j<ninputs; j++) {
	    if (matchval == INPUT_PARAM[j]) 
	      CHANGE_STAT[j] += s;
	  }
	}
      }
}
/*****************
 changestat: d_nodematch_sum
*****************/
WtC_CHANGESTAT_FN(c_nodematch_sum) { 
  Vertex ninputs = N_INPUT_PARAMS - N_NODES;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS();
      double matchval = INPUT_PARAM[tail+ninputs-1];
      if (matchval == INPUT_PARAM[head+ninputs-1]) { /* We have a match! */
	double s = weight - GETWT(tail,head);
	if (ninputs==0) {/* diff=F in network statistic specification */
	  CHANGE_STAT[0] += s;
	} else { /* diff=T */
	  for (unsigned int j=0; j<ninputs; j++) {
	    if (matchval == INPUT_PARAM[j]) 
	      CHANGE_STAT[j] += s;
	  }
	}
      }
}
/*****************
 changestat: d_nodemix_nonzero
 Update mixing matrix, non-bipartite networks only 
 (but see also d_mix_nonzero)
*****************/
WtC_CHANGESTAT_FN(c_nodemix_nonzero) {
  int ninputs = N_INPUT_PARAMS - N_NODES, ninputs2 = ninputs/2;
  double rtype, ctype, tmp;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS();
      double change = (weight!=0) - (GETWT(tail,head)!=0);
      /*Find the node covariate values (types) for the tail and head*/
      rtype=INPUT_PARAM[tail+ninputs-1];
      ctype=INPUT_PARAM[head+ninputs-1];
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
}
/*****************
 changestat: d_nodemix_sum
 Update mixing matrix, non-bipartite networks only 
 (but see also d_mix_sum)
*****************/
WtC_CHANGESTAT_FN(c_nodemix_sum) {
  int ninputs = N_INPUT_PARAMS - N_NODES, ninputs2 = ninputs/2;
  double rtype, ctype, tmp;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS();
      double change = weight - GETWT(tail,head);
      /*Find the node covariate values (types) for the tail and head*/
      rtype=INPUT_PARAM[tail+ninputs-1];
      ctype=INPUT_PARAM[head+ninputs-1];
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
}
/*****************
 stat: node o[ut] covar[iance] 
*****************/
WtC_CHANGESTAT_FN(c_nodeocovar){
  unsigned int transcode = INPUT_ATTRIB[0], center = INPUT_ATTRIB[1];
  double sum = center?*(double *)STORAGE : 0;
  
  ZERO_ALL_CHANGESTATS();
      double diff = TRANSFORM_DYADVAL(weight,transcode)-TRANSFORM_DYADVAL(GETWT(tail,head),transcode);
      double new_sum = center? sum + diff : 0;
      diff /= (N_NODES-2);
      EXEC_THROUGH_OUTEDGES(tail, e, i, yti, {
	  if(i!=head)
	    CHANGE_STAT[0] += 2*diff*TRANSFORM_DYADVAL(yti,transcode);
	});
      if(center){
	CHANGE_STAT[0] += (sum*sum-new_sum*new_sum)/N_DYADS;
	sum = new_sum;
      }
}
WtI_CHANGESTAT_FN(i_nodeocovar){
  unsigned int transcode = INPUT_ATTRIB[0], center = INPUT_ATTRIB[1];
  if(center){
    ALLOC_STORAGE(1, double, sum);
    *sum = 0;
    EXEC_THROUGH_NET_EDGES(tail, e1, head, y, {
	*sum+=TRANSFORM_DYADVAL(y, transcode);
	    head=head; e1=e1; // Prevent a compiler warning.
      });
  }
}

WtU_CHANGESTAT_FN(u_nodeocovar){
  unsigned int transcode = INPUT_ATTRIB[0], center = INPUT_ATTRIB[1];
  if(center){
    GET_STORAGE(double, sum);
    
    if(tail) *sum += TRANSFORM_DYADVAL(weight,transcode)-TRANSFORM_DYADVAL(GETWT(tail,head),transcode);
  }
}

/*****************
 stat: nodeofactor (nonzero)
*****************/
WtC_CHANGESTAT_FN(c_nodeofactor_nonzero){ 
  double s = (weight!=0) - (GETWT(tail,head)!=0);
  int tailpos = INPUT_ATTRIB[tail-1];
  if (tailpos != -1) CHANGE_STAT[tailpos] += s;
}
/*****************
 stat: nodeofactor (sum)
*****************/
WtC_CHANGESTAT_FN(c_nodeofactor_sum){ 
  double s = weight - GETWT(tail,head);
  int tailpos = INPUT_ATTRIB[tail-1];
  if (tailpos != -1) CHANGE_STAT[tailpos] += s;
}

/*****************
 stat: node sq[uare]r[oo]t covar[iance] "centered"
*****************/

WtC_CHANGESTAT_FN(c_nodesqrtcovar_centered){
  // Compute sum(sqrt(y)) (directed) or twice that (undirected). We can update it for each toggle.
  double ssq = *(double *)STORAGE;

  ZERO_ALL_CHANGESTATS();
      double change = 0;
      double sqrtdiff = sqrt(weight)-sqrt(GETWT(tail,head));
      double new_ssq = ssq + sqrtdiff*(DIRECTED? 1 : 2);
      
      EXEC_THROUGH_EDGES(tail, e, j, ytj, {
	  if(j!=head) change += sqrtdiff*sqrt(ytj);
	});
      
      EXEC_THROUGH_EDGES(head, e, i, yih, {
	  if(i!=tail) change += sqrtdiff*sqrt(yih);
	});

      CHANGE_STAT[0] += change/(N_NODES-2);
      CHANGE_STAT[0] -= (new_ssq*new_ssq-ssq*ssq) / (N_NODES*(N_NODES-1)) / 2;

      ssq = new_ssq;
}
WtU_CHANGESTAT_FN(u_nodesqrtcovar_centered){
  double *ssq;
  if(!STORAGE){
    STORAGE = Calloc(1, double);
    ssq = (double *)STORAGE;
    *ssq = 0;
    EXEC_THROUGH_NET_EDGES(i, e, j, yij, {
	*ssq += sqrt(yij); j=j; e=e; /* j=j and e=e are silly, just to prevent compiler warnings */
      });

  if(!DIRECTED) *ssq *= 2;
  }else ssq = (double *)STORAGE; 

  if(tail){
    double sqrtdiff = sqrt(weight)-sqrt(GETWT(tail,head));
    *ssq += sqrtdiff*(DIRECTED? 1 : 2);
  }
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
WtC_CHANGESTAT_FN(c_nodesqrtcovar){
  ZERO_ALL_CHANGESTATS();
      double change = 0;
      double sqrtdiff = sqrt(weight)-sqrt(GETWT(tail,head));
      
      EXEC_THROUGH_EDGES(tail, e, j, ytj, {
	  if(j!=head) change += sqrtdiff*sqrt(ytj);
	});
      
      EXEC_THROUGH_EDGES(head, e, i, yih, {
	  if(i!=tail) change += sqrtdiff*sqrt(yih);
	});

      CHANGE_STAT[0] += change/(N_NODES-2);
}
/*****************
 stat: nonzero
*****************/
WtC_CHANGESTAT_FN(c_nonzero){
  ZERO_ALL_CHANGESTATS();
        CHANGE_STAT[0] += (weight!=0) - (GETWT(tail,head)!=0);
}

/********************  changestats:   S    ***********/

/*****************
 stat: smallerthan
*****************/
WtC_CHANGESTAT_FN(c_smallerthan){
  ZERO_ALL_CHANGESTATS();
  double oldwt = GETWT(tail,head);
  for(unsigned int i=0; i<N_CHANGE_STATS; i++)
    CHANGE_STAT[i] += (weight<INPUT_ATTRIB[i]) - (oldwt<INPUT_ATTRIB[i]);
}

/*****************
 stat: sum
*****************/
WtC_CHANGESTAT_FN(c_sum){
  ZERO_ALL_CHANGESTATS();
      CHANGE_STAT[0] += weight-GETWT(tail,head);
}

/*****************
 stat: sum(with power)
*****************/
WtC_CHANGESTAT_FN(c_sum_pow){
  double p = INPUT_ATTRIB[0];
  
  ZERO_ALL_CHANGESTATS();
      CHANGE_STAT[0] += pow(weight,p)-pow(GETWT(tail,head),p);
}
/********************  changestats:   T    ***********/

/*****************
 stat: transitiveweights
*****************/

WtC_CHANGESTAT_FN(c_transitiveweights){ 
  CHANGE_STAT[0]=0;
  unsigned int path = INPUT_ATTRIB[0], combine = INPUT_ATTRIB[1], compare =  INPUT_ATTRIB[2];
  
  ZERO_ALL_CHANGESTATS();
      // (tail,head) as the focus dyad
      // This means that the combined 2-path value doesn't change.
      double two_paths = 0;
      EXEC_THROUGH_INEDGES(head, e1, k, ykh, {
	  double two_path;

	  switch(path){
	  case 1: two_path = fmin(GETWT(tail,k), ykh); break; // min
	  case 2: two_path = sqrt(GETWT(tail,k) * ykh); break; // geomean
	  default: // never reached, but prevents a warning
	    two_path = 0;    
	  }
	  
	  switch(combine){
	  case 1: two_paths = fmax(two_paths, two_path); break; // max
	  case 2: two_paths += two_path; break; // sum
	  }

	});

      switch(compare){
      case 1: CHANGE_STAT[0] += fmin(two_paths, weight) - fmin(two_paths, GETWT(tail,head)); break; // min
      case 2: CHANGE_STAT[0] += sqrt(two_paths * weight) - sqrt(two_paths * GETWT(tail,head)); break; // geomean
      }

      // (tail,head) as the first link of the 2-path
      // This means that only the combined 2-path value may change.
      EXEC_THROUGH_OUTEDGES(tail, e1, j, ytj, {
	  if(j==head) continue;
	  
	  double old_two_paths = 0;
	  double new_two_paths = 0;

	  EXEC_THROUGH_INEDGES(j, e2, k, ykj, {
	      double old_ytk = (k==head) ? GETWT(tail,head) : GETWT(tail, k);
	      double new_ytk = (k==head) ? weight : old_ytk; 
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

      // (tail,head) as the second link of the 2-path
      // This means that only the combined 2-path value may change.
      EXEC_THROUGH_INEDGES(head, e1, i, yih, {
	  if(i==tail) continue;
	  
	  double old_two_paths = 0;
	  double new_two_paths = 0;

	  EXEC_THROUGH_OUTEDGES(i, e2, k, yik, {
	      double old_ykh = (k==tail) ? GETWT(tail,head) : GETWT(k,head);
	      double new_ykh = (k==tail) ? weight : old_ykh; 
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

/*****************
 changestat: d_mixmat
 General mixing matrix (mm) implementation.
*****************/
WtC_CHANGESTAT_FN(c_mixmat_sum){
  unsigned int symm = ((int)INPUT_PARAM[0]) & 1;
  unsigned int marg = ((int)INPUT_PARAM[0]) & 2;
  double *tx = INPUT_PARAM;
  double *hx = BIPARTITE? INPUT_PARAM : INPUT_PARAM + N_NODES;
  double *cells = BIPARTITE? INPUT_PARAM + N_NODES + 1: INPUT_PARAM + N_NODES*2 + 1;
  
  unsigned int diag = tx[tail]==tx[head] && hx[tail]==hx[head];
  for(unsigned int j=0; j<N_CHANGE_STATS; j++){
    unsigned int thmatch = tx[tail]==cells[j*2] && hx[head]==cells[j*2+1];
    unsigned int htmatch = tx[head]==cells[j*2] && hx[tail]==cells[j*2+1];
    
    int w = DIRECTED || BIPARTITE? thmatch :
      (symm ? thmatch||htmatch : thmatch+htmatch)*(symm && marg && diag?2:1);
      if(w) CHANGE_STAT[j] += w*(weight-GETWT(tail,head));
  }
}

WtC_CHANGESTAT_FN(c_mixmat_nonzero){
  unsigned int symm = ((int)INPUT_PARAM[0]) & 1;
  unsigned int marg = ((int)INPUT_PARAM[0]) & 2;
  double *tx = INPUT_PARAM;
  double *hx = BIPARTITE? INPUT_PARAM : INPUT_PARAM + N_NODES;
  double *cells = BIPARTITE? INPUT_PARAM + N_NODES + 1: INPUT_PARAM + N_NODES*2 + 1;
  
  unsigned int diag = tx[tail]==tx[head] && hx[tail]==hx[head];
  for(unsigned int j=0; j<N_CHANGE_STATS; j++){
    unsigned int thmatch = tx[tail]==cells[j*2] && hx[head]==cells[j*2+1];
    unsigned int htmatch = tx[head]==cells[j*2] && hx[tail]==cells[j*2+1];
    
    int w = DIRECTED || BIPARTITE? thmatch :
      (symm ? thmatch||htmatch : thmatch+htmatch)*(symm && marg && diag?2:1);
    if(w) CHANGE_STAT[j] += w*((weight!=0)-(GETWT(tail,head)!=0));
  }
}
