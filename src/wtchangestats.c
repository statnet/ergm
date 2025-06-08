/*  File src/wtchangestats.c in package ergm, part of the Statnet suite of
 *  packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#include "wtchangestats.h"

/********************  changestats:   A    ***********/

/*****************
 stat: atleast
*****************/
WtC_CHANGESTAT_FN(c_atleast){
  for(unsigned int i=0; i<N_CHANGE_STATS; i++)
    CHANGE_STAT[i] += (weight>=INPUT_ATTRIB[i]) - (edgestate>=INPUT_ATTRIB[i]);
}

/*****************
 stat: atmost
*****************/
WtC_CHANGESTAT_FN(c_atmost){
  for(unsigned int i=0; i<N_CHANGE_STATS; i++)
    CHANGE_STAT[i] += (weight<=INPUT_ATTRIB[i]) - (edgestate<=INPUT_ATTRIB[i]);
}

/********************  changestats:   C    ***********/

/*****************
 stat: cyclicalweights
*****************/

WtC_CHANGESTAT_FN(c_cyclicalweights){ 
  unsigned int path = INPUT_ATTRIB[0], combine = INPUT_ATTRIB[1], compare =  INPUT_ATTRIB[2];
  CHANGE_STAT[0]=0;

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
      case 1: CHANGE_STAT[0] += fmin(two_paths, weight) - fmin(two_paths, edgestate); break; // min
      case 2: CHANGE_STAT[0] += sqrt(two_paths * weight) - sqrt(two_paths * edgestate); break; // geomean
      }

      // (tail,head) as the first link in the 2-path
      // This means that only the strongest 2-path may change.
      EXEC_THROUGH_INEDGES(tail, e1, j, yjt, {
	  if(j==head) continue;
	  
	  double old_two_paths = 0;
	  double new_two_paths = 0;

	  EXEC_THROUGH_INEDGES(j, e2, k, ykj, {
	      double old_ytk = (k==head) ? edgestate : GETWT(tail, k);
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
	      double old_ykh = (k==tail) ? edgestate : GETWT(k,head);
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
 stat: greaterthan
*****************/
WtC_CHANGESTAT_FN(c_greaterthan){
  for(unsigned int i=0; i<N_CHANGE_STATS; i++)
    CHANGE_STAT[i] += (weight>INPUT_ATTRIB[i]) - (edgestate>INPUT_ATTRIB[i]);
}

/********************  changestats:   I    ***********/

/*****************
 stat: ininterval
*****************/
WtC_CHANGESTAT_FN(c_ininterval){
      CHANGE_STAT[0] += ((INPUT_ATTRIB[2] ? weight>INPUT_ATTRIB[0] : weight>=INPUT_ATTRIB[0]) && (INPUT_ATTRIB[3] ? weight<INPUT_ATTRIB[1] : weight<=INPUT_ATTRIB[1])) - ((INPUT_ATTRIB[2] ? edgestate>INPUT_ATTRIB[0] : edgestate>=INPUT_ATTRIB[0]) && (INPUT_ATTRIB[3] ? edgestate<INPUT_ATTRIB[1] : edgestate<=INPUT_ATTRIB[1]));
}
/********************  changestats:   L    ***********/

/********************  changestats:   M    ***********/

/*****************
 stat: mutual (product a.k.a. covariance)
*****************/
WtC_CHANGESTAT_FN(c_mutual_wt_product){
      double htweight = GETWT(head,tail);
      CHANGE_STAT[0] += (weight*htweight) - (edgestate*htweight);
}
/*****************
 stat: mutual (geometric mean)
*****************/
WtC_CHANGESTAT_FN(c_mutual_wt_geom_mean){
      double htweight = GETWT(head,tail);
      CHANGE_STAT[0] += sqrt(weight*htweight) - sqrt(edgestate*htweight);
}
/*****************
 stat: mutual (minimum)
*****************/
WtC_CHANGESTAT_FN(c_mutual_wt_min){
      double htweight = GETWT(head,tail);
      CHANGE_STAT[0] += fmin(weight,htweight) - fmin(edgestate,htweight);
}

/*****************
 stat: mutual (-|yij-yji|)
*****************/
WtC_CHANGESTAT_FN(c_mutual_wt_nabsdiff){
      double htweight = GETWT(head,tail);
      CHANGE_STAT[0] -= fabs(weight-htweight) - fabs(edgestate-htweight);
}

/********************  changestats:   N    ***********/

/*****************
 stat: node covar[iance] 
*****************/
WtC_CHANGESTAT_FN(c_nodecovar){
  unsigned int transcode = INPUT_ATTRIB[0], center = INPUT_ATTRIB[1];
  double sum = center?*(double *)STORAGE : 0;
  
      double diff = TRANSFORM_DYADVAL(weight,transcode)-TRANSFORM_DYADVAL(edgestate,transcode);
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
    if(tail) *sum += TRANSFORM_DYADVAL(weight,transcode)-TRANSFORM_DYADVAL(edgestate,transcode);
  }
}

/*****************
 stat: node i[n] covar 
*****************/
WtC_CHANGESTAT_FN(c_nodeicovar){
  unsigned int transcode = INPUT_ATTRIB[0], center = INPUT_ATTRIB[1];
  double sum = center?*(double *)STORAGE : 0;
  
      double diff = TRANSFORM_DYADVAL(weight,transcode)-TRANSFORM_DYADVAL(edgestate,transcode);
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
    
    if(tail) *sum += TRANSFORM_DYADVAL(weight,transcode)-TRANSFORM_DYADVAL(edgestate,transcode);
  }
}

/*****************
 stat: node o[ut] covar[iance] 
*****************/
WtC_CHANGESTAT_FN(c_nodeocovar){
  unsigned int transcode = INPUT_ATTRIB[0], center = INPUT_ATTRIB[1];
  double sum = center?*(double *)STORAGE : 0;
  
      double diff = TRANSFORM_DYADVAL(weight,transcode)-TRANSFORM_DYADVAL(edgestate,transcode);
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
    
    if(tail) *sum += TRANSFORM_DYADVAL(weight,transcode)-TRANSFORM_DYADVAL(edgestate,transcode);
  }
}

/*****************
 stat: node sq[uare]r[oo]t covar[iance] "centered"
*****************/

WtC_CHANGESTAT_FN(c_nodesqrtcovar_centered){
  // Compute sum(sqrt(y)) (directed) or twice that (undirected). We can update it for each toggle.
  double ssq = *(double *)STORAGE;

      double change = 0;
      double sqrtdiff = sqrt(weight)-sqrt(edgestate);
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
    STORAGE = R_Calloc(1, double);
    ssq = (double *)STORAGE;
    *ssq = 0;
    EXEC_THROUGH_NET_EDGES(i, e, j, yij, {
	*ssq += sqrt(yij); j=j; e=e; /* j=j and e=e are silly, just to prevent compiler warnings */
      });

  if(!DIRECTED) *ssq *= 2;
  }else ssq = (double *)STORAGE; 

  if(tail){
    double sqrtdiff = sqrt(weight)-sqrt(edgestate);
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
      double change = 0;
      double sqrtdiff = sqrt(weight)-sqrt(edgestate);
      
      EXEC_THROUGH_EDGES(tail, e, j, ytj, {
	  if(j!=head) change += sqrtdiff*sqrt(ytj);
	});
      
      EXEC_THROUGH_EDGES(head, e, i, yih, {
	  if(i!=tail) change += sqrtdiff*sqrt(yih);
	});

      CHANGE_STAT[0] += change/(N_NODES-2);
}

/********************  changestats:   S    ***********/

/*****************
 stat: smallerthan
*****************/
WtC_CHANGESTAT_FN(c_smallerthan){
  for(unsigned int i=0; i<N_CHANGE_STATS; i++)
    CHANGE_STAT[i] += (weight<INPUT_ATTRIB[i]) - (edgestate<INPUT_ATTRIB[i]);
}

/*****************
 stat: sum(with power)
*****************/
WtC_CHANGESTAT_FN(c_sum_pow){
  double p = INPUT_ATTRIB[0];
  
      CHANGE_STAT[0] += pow(weight,p)-pow(edgestate,p);
}
/********************  changestats:   T    ***********/

/*****************
 stat: transitiveweights
*****************/

WtC_CHANGESTAT_FN(c_transitiveweights){ 
  CHANGE_STAT[0]=0;
  unsigned int path = INPUT_ATTRIB[0], combine = INPUT_ATTRIB[1], compare =  INPUT_ATTRIB[2];
  
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
      case 1: CHANGE_STAT[0] += fmin(two_paths, weight) - fmin(two_paths, edgestate); break; // min
      case 2: CHANGE_STAT[0] += sqrt(two_paths * weight) - sqrt(two_paths * edgestate); break; // geomean
      }

      // (tail,head) as the first link of the 2-path
      // This means that only the combined 2-path value may change.
      EXEC_THROUGH_OUTEDGES(tail, e1, j, ytj, {
	  if(j==head) continue;
	  
	  double old_two_paths = 0;
	  double new_two_paths = 0;

	  EXEC_THROUGH_INEDGES(j, e2, k, ykj, {
	      double old_ytk = (k==head) ? edgestate : GETWT(tail, k);
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
	      double old_ykh = (k==tail) ? edgestate : GETWT(k,head);
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
