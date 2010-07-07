#include "changestats.h"

/*****************
 changestat: d_degcor
*****************/
D_CHANGESTAT_FN(d_degcor) { 
  int i, echange;
  Vertex h, t, hdeg, tdeg, node3;
  Edge e;
  double sigma2;

  sigma2 = INPUT_PARAM[0];
// Rprintf("sigma2 %f\n",sigma2);
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i) {
    h = heads[i];
    t = tails[i];
    hdeg = OUT_DEG[h] + IN_DEG[h];
    tdeg = OUT_DEG[t] + IN_DEG[t];
    echange = IS_OUTEDGE(h, t) ? -1 : 1;
    if(echange==1){
     CHANGE_STAT[0] += (hdeg + 1.0)*(tdeg + 1.0);
     STEP_THROUGH_OUTEDGES(t, e, node3) { /* step through outedges of tail */
       CHANGE_STAT[0] += (OUT_DEG[node3] + IN_DEG[node3]);
     }
     STEP_THROUGH_INEDGES(t, e, node3) { /* step through inedges of tail */
       CHANGE_STAT[0] += (OUT_DEG[node3] + IN_DEG[node3]);
     }
     STEP_THROUGH_OUTEDGES(h, e, node3) { /* step through outedges of head */
       CHANGE_STAT[0] += (OUT_DEG[node3] + IN_DEG[node3]);
     }
     STEP_THROUGH_INEDGES(h, e, node3) { /* step through inedges of head */
       CHANGE_STAT[0] += (OUT_DEG[node3] + IN_DEG[node3]);
     }
    }else{
     CHANGE_STAT[0] -= (hdeg)*(tdeg);
     STEP_THROUGH_OUTEDGES(t, e, node3) { /* step through outedges of tail */
      if(node3!=h) CHANGE_STAT[0] -= (OUT_DEG[node3] + IN_DEG[node3]);
     }
     STEP_THROUGH_INEDGES(t, e, node3) { /* step through inedges of tail */
      if(node3!=h) CHANGE_STAT[0] -= (OUT_DEG[node3] + IN_DEG[node3]);
     }
     STEP_THROUGH_OUTEDGES(h, e, node3) { /* step through outedges of head */
      if(node3!=t) CHANGE_STAT[0] -= (OUT_DEG[node3] + IN_DEG[node3]);
     }
     STEP_THROUGH_INEDGES(h, e, node3) { /* step through inedges of head */
      if(node3!=t) CHANGE_STAT[0] -= (OUT_DEG[node3] + IN_DEG[node3]);
     }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
  CHANGE_STAT[0] *= (2.0/sigma2);
}
S_CHANGESTAT_FN(s_degcor) { 
  Vertex h, t, hdeg, tdeg;
  Edge e;
  double mu, mu2, sigma2, cross;

  mu = 0.0;
  mu2 = 0.0;
  cross = 0.0;
  for(h=1; h <= N_NODES; h++) {
   STEP_THROUGH_OUTEDGES(h, e, t) { /* step through outedges of head */
    hdeg = OUT_DEG[h] + IN_DEG[h];
    tdeg = OUT_DEG[t] + IN_DEG[t];
  // Rprintf("h %d t %d hdeg %d tdeg %d\n",h,t,hdeg,tdeg);
    mu  += hdeg + tdeg;
    mu2 += hdeg*hdeg + tdeg*tdeg;
    cross += 2.0*hdeg*tdeg;
   }
  }
  mu = mu / (2.0*N_EDGES);
  sigma2 = mu2/(2.0*N_EDGES) -  mu*mu;
// Rprintf("mu %f mu2 %f cross %f sigma2 %f\n",mu,mu2,cross,sigma2);
  CHANGE_STAT[0] = (cross / (2.0*N_EDGES) -  mu*mu) / sigma2;
}

/*****************
 changestat: d_adegcor
*****************/
D_CHANGESTAT_FN(d_adegcor) { 
  int i;
  double current;

  (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */
  current = mtp->dstats[0];
  FOR_EACH_TOGGLE(i) { TOGGLE(heads[i],tails[i]); }
  (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */
//  CHANGE_STAT[0] = mtp->dstats[0] - current;
//   Rprintf("c %f p %f",current,mtp->dstats[0]);
  mtp->dstats[0] -= current;
//   Rprintf(" p-c %f\n",mtp->dstats[0]);
  FOR_EACH_TOGGLE(i) { TOGGLE(heads[i],tails[i]); }
}
S_CHANGESTAT_FN(s_adegcor) { 
  Vertex h, t, hdeg, tdeg;
  Edge e;
  double mu, mu2, sigma2, cross;

  mu = 0.0;
  mu2 = 0.0;
  cross = 0.0;
  for(h=1; h <= N_NODES; h++) {
   STEP_THROUGH_OUTEDGES(h, e, t) { /* step through outedges of head */
    hdeg = OUT_DEG[h] + IN_DEG[h];
    tdeg = OUT_DEG[t] + IN_DEG[t];
    mu  += (double)(hdeg + tdeg);
    mu2 += (double)(hdeg*hdeg + tdeg*tdeg);
    cross += 2.0*hdeg*tdeg;
   }
  }
  mu = mu / (2.0*N_EDGES);
  sigma2 = mu2/(2.0*N_EDGES) -  mu*mu;
  CHANGE_STAT[0] = (cross / (2.0*N_EDGES) -  mu*mu) / sigma2;
}
D_CHANGESTAT_FN(d_rdegcor) { 
  int i;
  double current;

  (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */
  current = mtp->dstats[0];
  FOR_EACH_TOGGLE(i) { TOGGLE(heads[i],tails[i]); }
  (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */
//  CHANGE_STAT[0] = mtp->dstats[0] - current;
//   Rprintf("c %f p %f",current,mtp->dstats[0]);
  mtp->dstats[0] -= current;
//   Rprintf(" p-c %f\n",mtp->dstats[0]);
  FOR_EACH_TOGGLE(i) { TOGGLE(heads[i],tails[i]); }
}
S_CHANGESTAT_FN(s_rdegcor) { 
  Vertex h, t, hdeg, tdeg;
  Edge e;
  double mu, mu2, sigma2, cross;
  Vertex hrank, trank;
  Vertex *ndeg=malloc(sizeof(Vertex)*(N_NODES+1));

  for(h=0; h <= N_NODES; h++) { ndeg[h]=0; }
  for(h=0; h < N_NODES; h++) {
   STEP_THROUGH_OUTEDGES(h, e, t) { /* step through outedges of head */
    hdeg = OUT_DEG[h] + IN_DEG[h];
    tdeg = OUT_DEG[t] + IN_DEG[t];
    ndeg[hdeg+1]++;
    ndeg[tdeg+1]++;
   }
  }
for(h=1; h <= N_NODES; h++) {
    ndeg[h] += ndeg[h-1];
}
// Rprintf("h  %d hdeg[h] %d \n",h,ndeg[h]);}

  mu = 0.0;
  mu2 = 0.0;
  cross = 0.0;
  for(h=1; h <= N_NODES; h++) {
   STEP_THROUGH_OUTEDGES(h, e, t) { /* step through outedges of head */
    hdeg = OUT_DEG[h] + IN_DEG[h];
    tdeg = OUT_DEG[t] + IN_DEG[t];
    hrank = (ndeg[hdeg+1]+ndeg[hdeg+2]+1)*0.5;
    trank = (ndeg[tdeg+1]+ndeg[tdeg+2]+1)*0.5;
    mu  += (double)(hrank + trank);
    mu2 += (double)(hrank*hrank + trank*trank);
    cross += 2.0*hrank*trank;
   }
  }
  mu = mu / (2.0*N_EDGES);
  sigma2 = mu2/(2.0*N_EDGES) -  mu*mu;
  CHANGE_STAT[0] = (cross / (2.0*N_EDGES) -  mu*mu) / sigma2;
  free(ndeg);
}
