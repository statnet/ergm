#include "changestats.h"

/*****************
 changestat: d_degcor
*****************/
D_CHANGESTAT_FN(d_degcor) { 
  int i, echange;
  Vertex tail, head, taildeg, headdeg, node3;
  Edge e;
  double sigma2;

  sigma2 = INPUT_PARAM[0];
// Rprintf("sigma2 %f\n",sigma2);
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i) {
    tail = TAIL(i);
    head = HEAD(i);
    taildeg = OUT_DEG[tail] + IN_DEG[tail];
    headdeg = OUT_DEG[head] + IN_DEG[head];
    echange = IS_OUTEDGE(tail, head) ? -1 : 1;
    if(echange==1){
     CHANGE_STAT[0] += (taildeg + 1.0)*(headdeg + 1.0);
     STEP_THROUGH_OUTEDGES(head, e, node3) { /* step through outedges of head */
       CHANGE_STAT[0] += (OUT_DEG[node3] + IN_DEG[node3]);
     }
     STEP_THROUGH_INEDGES(head, e, node3) { /* step through inedges of head */
       CHANGE_STAT[0] += (OUT_DEG[node3] + IN_DEG[node3]);
     }
     STEP_THROUGH_OUTEDGES(tail, e, node3) { /* step through outedges of tail */
       CHANGE_STAT[0] += (OUT_DEG[node3] + IN_DEG[node3]);
     }
     STEP_THROUGH_INEDGES(tail, e, node3) { /* step through inedges of tail */
       CHANGE_STAT[0] += (OUT_DEG[node3] + IN_DEG[node3]);
     }
    }else{
     CHANGE_STAT[0] -= (taildeg)*(headdeg);
     STEP_THROUGH_OUTEDGES(head, e, node3) { /* step through outedges of head */
      if(node3!=tail) CHANGE_STAT[0] -= (OUT_DEG[node3] + IN_DEG[node3]);
     }
     STEP_THROUGH_INEDGES(head, e, node3) { /* step through inedges of head */
      if(node3!=tail) CHANGE_STAT[0] -= (OUT_DEG[node3] + IN_DEG[node3]);
     }
     STEP_THROUGH_OUTEDGES(tail, e, node3) { /* step through outedges of tail */
      if(node3!=head) CHANGE_STAT[0] -= (OUT_DEG[node3] + IN_DEG[node3]);
     }
     STEP_THROUGH_INEDGES(tail, e, node3) { /* step through inedges of tail */
      if(node3!=head) CHANGE_STAT[0] -= (OUT_DEG[node3] + IN_DEG[node3]);
     }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
  CHANGE_STAT[0] *= (2.0/sigma2);
}
S_CHANGESTAT_FN(s_degcor) { 
  Vertex tail, head, taildeg, headdeg;
  Edge e;
  double mu, mu2, sigma2, cross;

  mu = 0.0;
  mu2 = 0.0;
  cross = 0.0;
  for(tail=1; tail <= N_NODES; tail++) {
   STEP_THROUGH_OUTEDGES(tail, e, head) { /* step through outedges of tail */
    taildeg = OUT_DEG[tail] + IN_DEG[tail];
    headdeg = OUT_DEG[head] + IN_DEG[head];
  // Rprintf("tail %d head %d taildeg %d headdeg %d\n",tail,head,taildeg,headdeg);
    mu  += taildeg + headdeg;
    mu2 += taildeg*taildeg + headdeg*headdeg;
    cross += 2.0*taildeg*headdeg;
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
  FOR_EACH_TOGGLE(i) { TOGGLE(TAIL(i),HEAD(i)); }
  (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */
//  CHANGE_STAT[0] = mtp->dstats[0] - current;
//   Rprintf("c %f p %f",current,mtp->dstats[0]);
  mtp->dstats[0] -= current;
//   Rprintf(" p-c %f\n",mtp->dstats[0]);
FOR_EACH_TOGGLE(i) { TOGGLE(TAIL(i), HEAD(i)); }
}
S_CHANGESTAT_FN(s_adegcor) { 
  Vertex tail, head, taildeg, headdeg;
  Edge e;
  double mu, mu2, sigma2, cross;

  mu = 0.0;
  mu2 = 0.0;
  cross = 0.0;
  for(tail=1; tail <= N_NODES; tail++) {
   STEP_THROUGH_OUTEDGES(tail, e, head) { /* step through outedges of tail */
    taildeg = OUT_DEG[tail] + IN_DEG[tail];
    headdeg = OUT_DEG[head] + IN_DEG[head];
    mu  += (double)(taildeg + headdeg);
    mu2 += (double)(taildeg*taildeg + headdeg*headdeg);
    cross += 2.0*taildeg*headdeg;
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
  FOR_EACH_TOGGLE(i) { TOGGLE(TAIL(i), HEAD(i)); }
  (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */
//  CHANGE_STAT[0] = mtp->dstats[0] - current;
//   Rprintf("c %f p %f",current,mtp->dstats[0]);
  mtp->dstats[0] -= current;
//   Rprintf(" p-c %f\n",mtp->dstats[0]);
  FOR_EACH_TOGGLE(i) { TOGGLE(TAIL(i), HEAD(i)); }
}
S_CHANGESTAT_FN(s_rdegcor) { 
  Vertex tail, head, taildeg, headdeg;
  Edge e;
  double mu, mu2, sigma2, cross;
  Vertex tailrank, headrank;
  Vertex *ndeg=malloc(sizeof(Vertex)*(N_NODES+1));

  for(tail=0; tail <= N_NODES; tail++) { ndeg[tail]=0; }
  for(tail=0; tail < N_NODES; tail++) {
   STEP_THROUGH_OUTEDGES(tail, e, head) { /* step through outedges of tail */
    taildeg = OUT_DEG[tail] + IN_DEG[tail];
    headdeg = OUT_DEG[head] + IN_DEG[head];
    ndeg[taildeg+1]++;
    ndeg[headdeg+1]++;
   }
  }
for(tail=1; tail <= N_NODES; tail++) {
    ndeg[tail] += ndeg[tail-1];
}
// Rprintf("tail  %d taildeg[tail] %d \n",tail,ndeg[tail]);}

  mu = 0.0;
  mu2 = 0.0;
  cross = 0.0;
  for(tail=1; tail <= N_NODES; tail++) {
   STEP_THROUGH_OUTEDGES(tail, e, head) { /* step through outedges of tail */
    taildeg = OUT_DEG[tail] + IN_DEG[tail];
    headdeg = OUT_DEG[head] + IN_DEG[head];
    tailrank = (ndeg[taildeg+1]+ndeg[taildeg+2]+1)*0.5;
    headrank = (ndeg[headdeg+1]+ndeg[headdeg+2]+1)*0.5;
    mu  += (double)(tailrank + headrank);
    mu2 += (double)(tailrank*tailrank + headrank*headrank);
    cross += 2.0*tailrank*headrank;
   }
  }
  mu = mu / (2.0*N_EDGES);
  sigma2 = mu2/(2.0*N_EDGES) -  mu*mu;
  CHANGE_STAT[0] = (cross / (2.0*N_EDGES) -  mu*mu) / sigma2;
  free(ndeg);
}
S_CHANGESTAT_FN(s_pdegcor) { 
  Vertex tail, head, taildeg, headdeg;
  Edge e;
  double mu, mu2, mutail, mutail2, sigma2, sigmatail2, cross;

  mu = 0.0;
  mu2 = 0.0;
  mutail = 0.0;
  mutail2 = 0.0;
  cross = 0.0;
  for(tail=1; tail <= N_NODES; tail++) {
   STEP_THROUGH_OUTEDGES(tail, e, head) { /* step through outedges of tail */
    taildeg = OUT_DEG[tail];
    headdeg = IN_DEG[head];
    mu   += (double)(headdeg);
    mutail  += (double)(taildeg);
    mu2 += (double)(headdeg*headdeg);
    mutail2 += (double)(taildeg*taildeg);
    cross += taildeg*headdeg;
   }
  }
  mu = mu / (N_EDGES);
  mutail = mutail / (N_EDGES);
  sigma2 = mu2/(N_EDGES) -  mu*mu;
  sigmatail2 = mutail2/(N_EDGES) -  mutail*mutail;
  CHANGE_STAT[0] = (cross / (N_EDGES) -  mutail*mu) / sqrt(sigma2*sigmatail2);
}
D_CHANGESTAT_FN(d_pdegcor) { 
  int i;
  double current;

  (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */
  current = mtp->dstats[0];
  FOR_EACH_TOGGLE(i) { TOGGLE(TAIL(i), HEAD(i)); }
  (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */
  mtp->dstats[0] -= current;
  FOR_EACH_TOGGLE(i) { TOGGLE(TAIL(i), HEAD(i)); }
}
