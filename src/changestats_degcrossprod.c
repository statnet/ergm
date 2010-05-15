#include "changestats.h"

/*****************
 changestat: d_degcrossprod
*****************/
D_CHANGESTAT_FN(d_degcrossprod) { 
  int i, echange;
  Vertex h, t, hdeg, tdeg, node3;
  Edge e;
  double nedges;

  nedges = INPUT_PARAM[0];

  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i) {
    h = heads[i];
    t = tails[i];
    echange = IS_OUTEDGE(h, t) ? -1 : 1;
    hdeg = OUT_DEG[h] + IN_DEG[h];
    tdeg = OUT_DEG[t] + IN_DEG[t];
    if(echange==1){
     CHANGE_STAT[0] += (hdeg + 1)*(tdeg + 1);
     STEP_THROUGH_OUTEDGES(t, e, node3) { /* step through outedges of tail */
      if(node3!=h) CHANGE_STAT[0] += (OUT_DEG[node3] + IN_DEG[node3]);
     }
     STEP_THROUGH_INEDGES(t, e, node3) { /* step through inedges of tail */
      if(node3!=h) CHANGE_STAT[0] += (OUT_DEG[node3] + IN_DEG[node3]);
     }
     STEP_THROUGH_OUTEDGES(h, e, node3) { /* step through outedges of head */
      if(node3!=t) CHANGE_STAT[0] += (OUT_DEG[node3] + IN_DEG[node3]);
     }
     STEP_THROUGH_INEDGES(h, e, node3) { /* step through inedges of head */
      if(node3!=t) CHANGE_STAT[0] += (OUT_DEG[node3] + IN_DEG[node3]);
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
// Rprintf("N_EDGES %d nedges %f \n",N_EDGES, nedges);
  CHANGE_STAT[0] /= nedges;
}

/*****************
 changestat: d_degcor
*****************/
D_CHANGESTAT_FN(d_degcor) { 
  int i, echange;
  Vertex h, t, hdeg, tdeg, node3;
  Edge e;
  double sigma2;

  sigma2 = INPUT_PARAM[0];
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i) {
    h = heads[i];
    t = tails[i];
    echange = IS_OUTEDGE(h, t) ? -1 : 1;
    hdeg = OUT_DEG[h] + IN_DEG[h];
    tdeg = OUT_DEG[t] + IN_DEG[t];
    if(echange==1){
     CHANGE_STAT[0] += (hdeg + 1)*(tdeg + 1);
     STEP_THROUGH_OUTEDGES(t, e, node3) { /* step through outedges of tail */
      if(node3!=h) CHANGE_STAT[0] += (OUT_DEG[node3] + IN_DEG[node3]);
     }
     STEP_THROUGH_INEDGES(t, e, node3) { /* step through inedges of tail */
      if(node3!=h) CHANGE_STAT[0] += (OUT_DEG[node3] + IN_DEG[node3]);
     }
     STEP_THROUGH_OUTEDGES(h, e, node3) { /* step through outedges of head */
      if(node3!=t) CHANGE_STAT[0] += (OUT_DEG[node3] + IN_DEG[node3]);
     }
     STEP_THROUGH_INEDGES(h, e, node3) { /* step through inedges of head */
      if(node3!=t) CHANGE_STAT[0] += (OUT_DEG[node3] + IN_DEG[node3]);
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
  CHANGE_STAT[0] /= sigma2;
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
