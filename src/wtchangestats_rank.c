#include "wtchangestats_rank.h"

WtD_CHANGESTAT_FN(d_edgecov_rank){
  IF_1_EGO_SWAPS_2_ALTERS({
      Vertex v1=t;
      for (Vertex v2=1; v2 <= N_NODES; v2++){
	if(v2==v1) continue;
	double v12_old=GETWT(v1,v2);
	double v12_new=(v2!=h1 && v2!=h2) ? v12_old : ( v2==h1 ? weights[0] : weights[1] );
	for (Vertex v3=1; v3 <= N_NODES; v3++){
	  if(v3==v2 || v3==v1 || 
	     (h1!=v2 && h1!=v3 && h2!=v2 && h2!=v3)) continue;
	  double v13_old=GETWT(v1,v3);
	  double v13_new=(v3!=h1 && v3!=h2) ? v13_old : ( v3==h1 ? weights[0] : weights[1] );
	  double v123_covdiff=INPUT_PARAM[(v1-1)*N_NODES + (v2-1)] - INPUT_PARAM[(v1-1)*N_NODES + (v3-1)];
	  if(v12_old>v13_old)
	    CHANGE_STAT[0] -= v123_covdiff;
	  if(v12_new>v13_new)
	    CHANGE_STAT[0] += v123_covdiff;
	}
      }
    });
} 

WtS_CHANGESTAT_FN(s_edgecov_rank){
  CHANGE_STAT[0]=0;
  for (Vertex v1=1; v1 <= N_NODES; v1++){
    for (Vertex v2=1; v2 <= N_NODES; v2++){
      if(v2==v1) continue;
      double v12=GETWT(v1,v2);
      for (Vertex v3=1; v3 <= N_NODES; v3++){
	if(v3==v2 || v3==v1) continue;
	if(v12>GETWT(v1,v3))
	  CHANGE_STAT[0] += INPUT_PARAM[(v1-1)*N_NODES + (v2-1)] - INPUT_PARAM[(v1-1)*N_NODES + (v3-1)];
      }
    }
  }
}


WtD_CHANGESTAT_FN(d_inconsistency_rank){
  IF_1_EGO_SWAPS_2_ALTERS({
      Vertex v1=t;
      for(Vertex v2=1; v2 <= N_NODES; v2++){
	if(v2==v1) continue;

	// For some reason completely beyond me, CPP complains if I
	// declare more than one variable in a line while inside a
	// macro.
	double v12_old = GETWT(v1,v2);
	double v12_ref = INPUT_PARAM[(v1-1)*N_NODES+(v2-1)];
	double v12_new = (v2!=h1 && v2!=h2) ? v12_old : ( v2==h1 ? weights[0] : weights[1] );
	for(Vertex v3=1; v3 <= N_NODES; v3++){
	  if(v3==v2 || v3==v1 || 
	     (h1!=v2 && h1!=v3 && h2!=v2 && h2!=v3)) continue;
	  double v13_old=GETWT(v1,v3);
	  double v13_ref=INPUT_PARAM[(v1-1)*N_NODES+(v3-1)];
	  double v13_new=(v3!=h1 && v3!=h2) ? v13_old : ( v3==h1 ? weights[0] : weights[1] );
	  if((v12_old>v13_old)!=(v12_ref>v13_ref)) CHANGE_STAT[0]--;
	  if((v12_new>v13_new)!=(v12_ref>v13_ref)) CHANGE_STAT[0]++;
	}
      }
    });
}

WtS_CHANGESTAT_FN(s_inconsistency_rank){ 
  CHANGE_STAT[0]=0;
  for(Vertex v1=1; v1 <= N_NODES; v1++){
    for(Vertex v2=1; v2 <= N_NODES; v2++){
      if(v2==v1) continue;
      double v12 = GETWT(v1,v2), v12_ref = INPUT_PARAM[(v1-1)*N_NODES+(v2-1)];
      for(Vertex v3=1; v3 <= N_NODES; v3++){
	if(v3==v2 || v3==v1) continue;
	unsigned int 
	  v123 = v12>GETWT(v1,v3),
	  v123_ref = v12_ref>INPUT_PARAM[(v1-1)*N_NODES+(v3-1)];
	if(v123!=v123_ref) CHANGE_STAT[0]++;
      }
    }
  }
}

WtD_CHANGESTAT_FN(d_inconsistency_cov_rank){
  IF_1_EGO_SWAPS_2_ALTERS({
      unsigned int cov_start = N_NODES*N_NODES;
      Vertex v1=t;
      for(Vertex v2=1; v2 <= N_NODES; v2++){
	if(v2==v1) continue;

	// For some reason completely beyond me, CPP complains if I
	// declare more than one variable in a line while inside a
	// macro.
	double v12_old = GETWT(v1,v2);
	double v12_ref = INPUT_PARAM[(v1-1)*N_NODES+(v2-1)];
	double v12_new = (v2!=h1 && v2!=h2) ? v12_old : ( v2==h1 ? weights[0] : weights[1] );
	for(Vertex v3=1; v3 <= N_NODES; v3++){
	  if(v3==v2 || v3==v1 || 
	     (h1!=v2 && h1!=v3 && h2!=v2 && h2!=v3)) continue;
	  double v13_old=GETWT(v1,v3);
	  double v13_ref=INPUT_PARAM[(v1-1)*N_NODES+(v3-1)];
	  double v13_new=(v3!=h1 && v3!=h2) ? v13_old : ( v3==h1 ? weights[0] : weights[1] );
	  if((v12_old>v13_old)!=(v12_ref>v13_ref)) CHANGE_STAT[0]-=INPUT_PARAM[cov_start + (v1-1)*N_NODES*N_NODES + (v2-1)*N_NODES + (v3-1)];
	  if((v12_new>v13_new)!=(v12_ref>v13_ref)) CHANGE_STAT[0]+=INPUT_PARAM[cov_start + (v1-1)*N_NODES*N_NODES + (v2-1)*N_NODES + (v3-1)];
	}
      }
    });
}

WtS_CHANGESTAT_FN(s_inconsistency_cov_rank){ 
  unsigned int cov_start = N_NODES*N_NODES;
  CHANGE_STAT[0]=0;
  for(Vertex v1=1; v1 <= N_NODES; v1++){
    for(Vertex v2=1; v2 <= N_NODES; v2++){
      if(v2==v1) continue;
      double v12 = GETWT(v1,v2), v12_ref = INPUT_PARAM[(v1-1)*N_NODES+(v2-1)];
      for(Vertex v3=1; v3 <= N_NODES; v3++){
	if(v3==v2 || v3==v1) continue;
	unsigned int 
	  v123 = v12>GETWT(v1,v3),
	  v123_ref = v12_ref>INPUT_PARAM[(v1-1)*N_NODES+(v3-1)];
	if(v123!=v123_ref) 
	  CHANGE_STAT[0]+=INPUT_PARAM[cov_start + (v1-1)*N_NODES*N_NODES + (v2-1)*N_NODES + (v3-1)];
      }
    }
  }
}

WtD_CHANGESTAT_FN(d_deference){
  IF_1_EGO_SWAPS_2_ALTERS({
      for(Vertex v1=1; v1 <= N_NODES; v1++){
	for(Vertex v3=1; v3 <= N_NODES; v3++){
	  if(v3==v1) continue;
	  double v31_old = GETWT(v3,v1);
	  double v13_old = GETWT(v1,v3);
	  double v31_new = (v3!=t || (v1!=h1 && v1!=h2)) ? v31_old :
	    ( v1==h1 ? weights[0] : weights[1] );
	  double v13_new = (v1!=t || (v3!=h1 && v3!=h2)) ? v13_old :
	    ( v3==h1 ? weights[0] : weights[1] );

	  for(Vertex v2=1; v2 <= N_NODES; v2++){
	    if(v2==v1 || v3==v2 || 
	       (t!=v1 && t!=v3 && 
		h1!=v1 && h1!=v2 && h1!=v3 &&
		h2!=v1 && h2!=v2 && h2!=v3)) continue;
	    double v32_old = GETWT(v3,v2);
	    double v12_old = GETWT(v1,v2);
	    double v32_new = (v3!=t || (v2!=h1 && v2!=h2)) ? v32_old :
	      ( v2==h1 ? weights[0] : weights[1] );
	    double v12_new = (v1!=t || (v2!=h1 && v2!=h2)) ? v12_old :
	      ( v2==h1 ? weights[0] : weights[1] );
	    if(v32_old>v31_old && v13_old>v12_old) 
	      CHANGE_STAT[0]--;
	    if(v32_new>v31_new && v13_new>v12_new) 
	      CHANGE_STAT[0]++;
	  }
	}
      }
    });
}

WtS_CHANGESTAT_FN(s_deference){ 
  CHANGE_STAT[0]=0;
  for(Vertex v1=1; v1 <= N_NODES; v1++){
    for(Vertex v3=1; v3 <= N_NODES; v3++){
      if(v3==v1) continue;
      double v31 = GETWT(v3,v1), v13 = GETWT(v1,v3);
      for(Vertex v2=1; v2 <= N_NODES; v2++){
	if(v2==v1 || v3==v2) continue;
	if(GETWT(v3,v2)>v31 && v13>GETWT(v1,v2)) 
	  CHANGE_STAT[0]++;
      }
    }
  }
}

WtD_FROM_S_FN(d_nodeicov_rank)

WtS_CHANGESTAT_FN(s_nodeicov_rank){ 
  CHANGE_STAT[0]=0;
  for (Vertex v1=1; v1 <= N_NODES; v1++){
    for (Vertex v2=1; v2 <= N_NODES; v2++){
      if(v2==v1) continue;
      double v12=GETWT(v1,v2);
      for (Vertex v3=1; v3 <= N_NODES; v3++){
	if(v3==v2 || v3==v1) continue;
	if(v12>GETWT(v1,v3)) 
	  CHANGE_STAT[0] += INPUT_PARAM[v2] - INPUT_PARAM[v3];
      }
    }
  }
}

WtD_CHANGESTAT_FN(d_nonconformity){
  IF_1_EGO_SWAPS_2_ALTERS({
      Vertex v1=t;
      
      for(Vertex v2=1; v2 <= N_NODES; v2++){
	if(v2==v1) continue;
	
	for(Vertex v3=1; v3 <= N_NODES; v3++){
	  if(v3==v2 || v3==v1) continue;
	  
	  double v13_old=GETWT(v1,v3);
	  double v23=GETWT(v2,v3);
	  double v13_new = (v3!=h1 && v3!=h2) ? v13_old :
	    ( v3==h1 ? weights[0] : weights[1] );
	  
	  for(Vertex v4=1; v4 <= N_NODES; v4++){
	    if(v4==v3 || v4==v2 || v4==v1 || 
	       (h1!=v3 && h1!=v4 && h2!=v3 && h2!=v4)) continue;
	    
	    double v14_old=GETWT(v1,v4);
	    double v24=GETWT(v2,v4);
	    double v14_new = (v4!=h1 && v4!=h2) ? v14_old :
	      ( v4==h1 ? weights[0] : weights[1] );
	    
	    if((v13_old>v14_old)!=(v23>v24)) CHANGE_STAT[0]-=2;
	    if((v13_new>v14_new)!=(v23>v24)) CHANGE_STAT[0]+=2;
	  }
	}
      }
    });
}

WtS_CHANGESTAT_FN(s_nonconformity){ 
  CHANGE_STAT[0]=0;
  for(Vertex v1=1; v1 <= N_NODES; v1++){
    for(Vertex v2=1; v2 < v1; v2++){
      for(Vertex v3=1; v3 <= N_NODES; v3++){
	if(v3==v2 || v3==v1) continue;
	double v13=GETWT(v1,v3), v23=GETWT(v2,v3);
	for(Vertex v4=1; v4 <= N_NODES; v4++){
	  if(v4==v3 || v4==v2 || v4==v1) continue;
	  unsigned int
	    v134 = v13>GETWT(v1,v4),
	    v234 = v23>GETWT(v2,v4);
	  if(v134!=v234) CHANGE_STAT[0]+=2;
	}
      }
    }
  }
}

WtD_CHANGESTAT_FN(d_local_nonconformity){
  IF_1_EGO_SWAPS_2_ALTERS({
      Vertex v1=t;
      for(Vertex v2=1; v2 <= N_NODES; v2++){
	if(v2==v1) continue;
	double v12_old=GETWT(v1,v2);
	double v12_new=(v1!=t || (v2!=h1 && v2!=h2)) ? v12_old :
	  ( v2==h1 ? weights[0] : weights[1] );
	for(Vertex v3=1; v3 <= N_NODES; v3++){
	  if(v3==v2 || v3==v1) continue;
	  double v13_old=GETWT(v1,v3);
	  double v13_new=(v1!=t || (v3!=h1 && v3!=h2)) ? v13_old :
	    ( v3==h1 ? weights[0] : weights[1] );
	  
	  if(v13_old<=v12_old && v13_new<=v12_new) continue;
	  
	  double v32_old=GETWT(v3,v2);
	  double v32_new=(v3!=t || (v2!=h1 && v2!=h2)) ? v32_old :
	    ( v2==h1 ? weights[0] : weights[1] );
	  
	  for(Vertex v4=1; v4 <= N_NODES; v4++){
	    if(v4==v3 || v4==v2 || v4==v1 ||
	       (h1!=v2 && h2!=v2 && h1!=v3 && h2!=v3 && h1!=v4 && h2!=v4)) continue;
	    double v14_old=GETWT(v1,v4);
	    double v14_new=(v1!=t || (v4!=h1 && v4!=h2)) ? v14_old :
	      ( v4==h1 ? weights[0] : weights[1] );
	    double v34_old=GETWT(v3,v4);
	    double v34_new=(v3!=t || (v4!=h1 && v4!=h2)) ? v34_old :
	      ( v4==h1 ? weights[0] : weights[1] );
	    
	    if(v13_old>v12_old && v14_old<=v12_old && v34_old>v32_old) CHANGE_STAT[0]--;
	    if(v13_new>v12_new && v14_new<=v12_new && v34_new>v32_new) CHANGE_STAT[0]++;
	  }
	}
      }
    
      for(Vertex v1=1; v1 <= N_NODES; v1++){
	for(Vertex v2=1; v2 <= N_NODES; v2++){
	  if(v2==v1) continue;
	  double v12_old=GETWT(v1,v2);
	  double v12_new=(v1!=t || (v2!=h1 && v2!=h2)) ? v12_old :
	    ( v2==h1 ? weights[0] : weights[1] );
	  Vertex v3=t;
	  if(v3==v2 || v3==v1) continue;
	  double v13_old=GETWT(v1,v3);
	  double v13_new=(v1!=t || (v3!=h1 && v3!=h2)) ? v13_old :
	    ( v3==h1 ? weights[0] : weights[1] );
	  
	  if(v13_old<=v12_old && v13_new<=v12_new) continue;
	  
	  double v32_old=GETWT(v3,v2);
	  double v32_new=(v3!=t || (v2!=h1 && v2!=h2)) ? v32_old :
	    ( v2==h1 ? weights[0] : weights[1] );
	  
	  for(Vertex v4=1; v4 <= N_NODES; v4++){
	    if(v4==v3 || v4==v2 || v4==v1 ||
	       (h1!=v2 && h2!=v2 && h1!=v3 && h2!=v3 && h1!=v4 && h2!=v4)) continue;
	    double v14_old=GETWT(v1,v4);
	    double v14_new=(v1!=t || (v4!=h1 && v4!=h2)) ? v14_old :
	      ( v4==h1 ? weights[0] : weights[1] );
	    double v34_old=GETWT(v3,v4);
	    double v34_new=(v3!=t || (v4!=h1 && v4!=h2)) ? v34_old :
		( v4==h1 ? weights[0] : weights[1] );
	    
	    if(v13_old>v12_old && v14_old<=v12_old && v34_old>v32_old) CHANGE_STAT[0]--;
	    if(v13_new>v12_new && v14_new<=v12_new && v34_new>v32_new) CHANGE_STAT[0]++;
	  }
	}
      }
    });
}

WtS_CHANGESTAT_FN(s_local_nonconformity){ 
  CHANGE_STAT[0]=0;
  for(Vertex v1=1; v1 <= N_NODES; v1++){
    for(Vertex v2=1; v2 <= N_NODES; v2++){
      if(v2==v1) continue;
      double v12=GETWT(v1,v2);
      for(Vertex v3=1; v3 <= N_NODES; v3++){
	if(v3==v2 || v3==v1) continue;
	double v13=GETWT(v1,v3);
	if(v13<=v12) continue;
	double v32=GETWT(v3,v2);
	for(Vertex v4=1; v4 <= N_NODES; v4++){
	  if(v4==v3 || v4==v2 || v4==v1) continue;
	  double v14=GETWT(v1,v4);
	  double v34=GETWT(v3,v4);
	  if(v14<=v12 && v34>v32) CHANGE_STAT[0]++;
	}
      }
    }
  }
}

WtD_FROM_S_FN(d_nonconformity_decay)

WtS_CHANGESTAT_FN(s_nonconformity_decay){ 
  CHANGE_STAT[0]=0;
  for(Vertex v1=1; v1 <= N_NODES; v1++){
    for(Vertex v2=1; v2 <= N_NODES; v2++){
      if(v2==v1) continue;
      double e = pow(INPUT_PARAM[1],INPUT_PARAM[0]-GETWT(v1,v2));
      for(Vertex v3=1; v3 <= N_NODES; v3++){
	if(v3==v2 || v3==v1) continue;
	double v13=GETWT(v1,v3), v23=GETWT(v2,v3);
	for(Vertex v4=1; v4 <= N_NODES; v4++){
	  if(v4==v3 || v4==v2 || v4==v1) continue;
	  unsigned int
	    v134 = v13>GETWT(v1,v4),
	    v234 = v23>GETWT(v2,v4);
	  if(v134!=v234) CHANGE_STAT[0]+=e;
	}
      }
    }
  }
}

WtD_FROM_S_FN(d_nonconformity_thresholds)

WtS_CHANGESTAT_FN(s_nonconformity_thresholds){ 
  ZERO_ALL_CHANGESTATS();
  for(Vertex v1=1; v1 <= N_NODES; v1++){
    for(Vertex v2=1; v2 <= N_NODES; v2++){
      if(v2==v1) continue;
      double v12 = GETWT(v1,v2);
      unsigned int i;
      for(i=0; i<N_CHANGE_STATS; i++) if(v12>=INPUT_PARAM[i]) break;
      if(i==N_CHANGE_STATS) continue;
      for(Vertex v3=1; v3 <= N_NODES; v3++){
	if(v3==v2 || v3==v1) continue;
	double v13=GETWT(v1,v3), v23=GETWT(v2,v3);
	for(Vertex v4=1; v4 <= N_NODES; v4++){
	  if(v4==v3 || v4==v2 || v4==v1) continue;
	  unsigned int
	    v134 = v13>GETWT(v1,v4),
	    v234 = v23>GETWT(v2,v4);
	  if(v134!=v234) CHANGE_STAT[i]+=1;
	}
      }
    }
  }
}

