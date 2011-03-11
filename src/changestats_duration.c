#include "changestats_duration.h"

/*****************
 void d_D_off
 This gives the change in mean off-duration for all non-existant edges
 (It does not yet work.)
*****************/
D_CHANGESTAT_FN(d_D_off) {
  int i;
  
  ZERO_ALL_CHANGESTATS(i); 
  FOR_EACH_TOGGLE(i) {
    TOGGLE_IF_MORE_TO_COME(i); 
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 void d_D_dyad
 This gives the change in mean duration for all dyads (where dyads without
 an edge are considered to have duration zero)
*****************/
D_CHANGESTAT_FN(d_D_dyad) {
  int edgeflag, i;
  Vertex tail, head;
  long int elapsed;

  ZERO_ALL_CHANGESTATS(i); 
  FOR_EACH_TOGGLE(i) {
    edgeflag = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i));
    if (edgeflag) {
      elapsed = ElapsedTime(tail,head,nwp);
      CHANGE_STAT[0] += (double) (N_EDGES - 1.0 - elapsed) / N_DYADS;
    } else {
      /* If we gain an edge, the change is always negative */
      CHANGE_STAT[0] += (double)(N_EDGES + 1.0) / N_DYADS;
    }
    TOGGLE_IF_MORE_TO_COME(i); 
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 void d_D_edge
 This gives the change in mean duration for existant edges
*****************/
D_CHANGESTAT_FN(d_D_edge) {
  int edgeflag, i;
  Vertex tail, head;
  Edge e=nwp->nedges;
  long int elapsed;
  double md = mean_duration(nwp);

  ZERO_ALL_CHANGESTATS(i); 
  FOR_EACH_TOGGLE(i) {
    edgeflag = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i));
    if (edgeflag) {
      elapsed = ElapsedTime(tail, head, nwp);
      /* Make sure to prevent division by zero */
      CHANGE_STAT[0] += e==1 ? -(double)elapsed : (md - elapsed)/(e-1.0);      
    } else {
      CHANGE_STAT[0] += (1.0-md)/(e+1.0);
    }
    TOGGLE_IF_MORE_TO_COME(i); 
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

double mean_duration(Network *nwp) 
{
  Vertex i, j;  
  Edge k, e=nwp->nedges;
  double sum=0.0;
    
  /* Repeated use of FindithEdge makes for inefficient code at the moment. */
  for (k=1; k <= e; k++) {
    FindithEdge(&i, &j, k, nwp);
    sum += ElapsedTime(i, j, nwp);
  }
  /* Protect against division by zero */
  return e==0 ? 0.0 : sum/(double)e;
}


/*****************
 void d_edges_ageinterval

 This is essentially the edges statistic, which only counts dyads with "age"
 (time steps spent in the current state) in the interval [inputparams0,inputparams1).
*****************/
D_CHANGESTAT_FN(d_edges_ageinterval){
  int edgeflag, i;
  Vertex tail, head;
  int from=mtp->inputparams[0], to=mtp->inputparams[1];
  
  ZERO_ALL_CHANGESTATS(i); 
  FOR_EACH_TOGGLE(i) {
    int age = ElapsedTime(tail=TAIL(i),head=HEAD(i),nwp);
    // Only count if the age is in [from,to). ( to=0 ==> to=Inf )
    if(from<=age && (to==0 || age<to)){
      edgeflag = IS_OUTEDGE(tail, head);
      CHANGE_STAT[0] += edgeflag ? - 1 : 1;
    }
    TOGGLE_IF_MORE_TO_COME(i); 
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

