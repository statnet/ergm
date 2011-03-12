#include "wtchangestats_rank.h"

WtD_FROM_S_FN(d_inconsistency_rank)

WtS_CHANGESTAT_FN(s_inconsistency_rank){ 
  CHANGE_STAT[0]=0;
  for(Vertex v1=1; v1 <= N_NODES; v1++){
    for(Vertex v2=1; v2 <= N_NODES; v2++){
      if(v2==v1) continue;
      for(Vertex v3=1; v3 <= N_NODES; v3++){
	if(v3==v2 || v3==v1) continue;
	unsigned int 
	  v123 = GETWT(v1,v2)>GETWT(v1,v3),
	  v123_ref = INPUT_PARAM[(v1-1)*N_NODES+(v2-1)]>INPUT_PARAM[(v1-1)*N_NODES+(v3-1)];
	if(v123!=v123_ref) CHANGE_STAT[0]++;
      }
    }
  }
}

WtD_FROM_S_FN(d_deference)

WtS_CHANGESTAT_FN(s_deference){ 
  CHANGE_STAT[0]=0;
  for(Vertex v1=1; v1 <= N_NODES; v1++){
    for(Vertex v2=1; v2 <= N_NODES; v2++){
      if(v2==v1) continue;
      for(Vertex v3=1; v3 <= N_NODES; v3++){
	if(v3==v2 || v3==v1) continue;
	if(GETWT(v3,v2)>GETWT(v3,v1) && GETWT(v1,v3)>GETWT(v1,v2)) 
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
      for (Vertex v3=1; v3 <= N_NODES; v3++){
	if(v3==v2 || v3==v1) continue;
	if(GETWT(v1,v2)>GETWT(v1,v3)) 
	  CHANGE_STAT[0] += INPUT_PARAM[v2] - INPUT_PARAM[v3];
      }
    }
  }
}

WtD_FROM_S_FN(d_nonconformity)

WtS_CHANGESTAT_FN(s_nonconformity){ 
  CHANGE_STAT[0]=0;
  for(Vertex v1=1; v1 <= N_NODES; v1++){
    for(Vertex v2=1; v2 < v1; v2++){
      for(Vertex v3=1; v3 <= N_NODES; v3++){
	if(v3==v2 || v3==v1) continue;
	for(Vertex v4=1; v4 <= N_NODES; v4++){
	  if(v4==v3 || v4==v2 || v4==v1) continue;
	  unsigned int
	    v134 = GETWT(v1,v3)>GETWT(v1,v4),
	    v234 = GETWT(v2,v3)>GETWT(v2,v4);
	  if(v134!=v234) CHANGE_STAT[0]+=2;
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
	for(Vertex v4=1; v4 <= N_NODES; v4++){
	  if(v4==v3 || v4==v2 || v4==v1) continue;
	  unsigned int
	    v134 = GETWT(v1,v3)>GETWT(v1,v4),
	    v234 = GETWT(v2,v3)>GETWT(v2,v4);
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
	for(Vertex v4=1; v4 <= N_NODES; v4++){
	  if(v4==v3 || v4==v2 || v4==v1) continue;
	  unsigned int
	    v134 = GETWT(v1,v3)>GETWT(v1,v4),
	    v234 = GETWT(v2,v3)>GETWT(v2,v4);
	  if(v134!=v234) CHANGE_STAT[i]+=1;
	}
      }
    }
  }
}

