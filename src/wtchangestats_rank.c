#include "wtchangestats_rank.h"

WtD_FROM_S_FN(d_deference)

WtS_CHANGESTAT_FN(s_deference){ 
  Vertex v1, v2, v3;
  
  CHANGE_STAT[0]=0;
  for (v1=1; v1 <= N_NODES; v1++){
    for (v2=1; v2 <= N_NODES; v2++){
      if(v2==v1) continue;
      for (v3=1; v3 <= N_NODES; v3++){
	if(v3==v2 || v3==v1) continue;
	if(GETWT(v3,v2)>GETWT(v3,v1) && GETWT(v1,v3)>GETWT(v1,v2)) 
	  CHANGE_STAT[0]++;
      }
    }
  }
}

WtD_FROM_S_FN(d_nodeicov_rank)

WtS_CHANGESTAT_FN(s_nodeicov_rank){ 
  Vertex v1, v2, v3;
  
  CHANGE_STAT[0]=0;
  for (v1=1; v1 <= N_NODES; v1++){
    for (v2=1; v2 <= N_NODES; v2++){
      if(v2==v1) continue;
      for (v3=1; v3 <= N_NODES; v3++){
	if(v3==v2 || v3==v1) continue;
	if(GETWT(v1,v2)>GETWT(v1,v3)) 
	  CHANGE_STAT[0] += INPUT_ATTRIB[v2] - INPUT_ATTRIB[v3];
      }
    }
  }
}

WtD_FROM_S_FN(d_nonconformity)

WtS_CHANGESTAT_FN(s_nonconformity){ 
  Vertex v1, v2, v3, v4;
  
  CHANGE_STAT[0]=0;
  for (v1=1; v1 <= N_NODES; v1++){
    for (v2=1; v2 <= N_NODES; v2++){
      if(v2==v1) continue;
      for (v3=1; v3 <= N_NODES; v3++){
	if(v3==v2 || v3==v1) continue;
	for (v4=1; v4 <= N_NODES; v4++){
	  if(v4==v3 || v4==v2 || v4==v1) continue;
	  unsigned int
	    v123 = GETWT(v1,v2)>GETWT(v1,v3)?1:0,
	    v423 = GETWT(v4,v2)>GETWT(v4,v3)?1:0;
	  CHANGE_STAT[0]+=v123*(1-v423)+(1-v123)*v423;
	}
      }
    }
  }
}
