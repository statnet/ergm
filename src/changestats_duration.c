/*
 *  File ergm/src/changestats_duration.c
 *  Part of the statnet package, http://statnet.org
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) in
 *    http://statnet.org/attribution
 *
 *  Copyright 2012 the statnet development team
 */
#include "changestats_duration.h"


/*****************
 void d_edges_ageinterval

 This is essentially the edges statistic, which only counts dyads with "age"
 (time steps spent in the current state) in the interval [inputparams0,inputparams1).

It changes whenever the clock advances, so it needs a t_??? function.
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

D_CHANGESTAT_FN(d_edges_ageinterval_mon){
  int edgeflag, i;
  Vertex tail, head;
  int from=mtp->inputparams[0], to=mtp->inputparams[1];
  
  ZERO_ALL_CHANGESTATS(i); 
  FOR_EACH_TOGGLE(i) {
    int age = ElapsedTime(tail=TAIL(i),head=HEAD(i),nwp);
    // Only count if the age is in [from,to). ( to=0 ==> to=Inf )
    if(from<=age && (to==0 || age<to)){
      edgeflag = IS_OUTEDGE(tail, head);
      CHANGE_STAT[0] += edgeflag ? - 1 : 0;
    }
    TOGGLE_IF_MORE_TO_COME(i); 
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

T_CHANGESTAT_FN(t_edges_ageinterval_mon){
  int from=mtp->inputparams[0], to=mtp->inputparams[1];
  CHANGE_STAT[0] = 0;
  for (Edge k=1; k <= N_EDGES; k++) {
    Vertex tail, head;
    FindithEdge(&tail, &head, k, nwp);
    int age = ElapsedTime(tail,head,nwp) + 1; // Every tie "starts out" at age 1.
    if(age == from) CHANGE_STAT[0]++; // The tie "ages" into the interval.
    if(to!=0 && age == to) CHANGE_STAT[0]--; // The tie "ages" out of the interval.
  }
}

S_CHANGESTAT_FN(s_edges_ageinterval_mon){
  int from=mtp->inputparams[0], to=mtp->inputparams[1];
  CHANGE_STAT[0] = 0;
  for (Edge k=1; k <= N_EDGES; k++) {
    Vertex tail, head;
    FindithEdge(&tail, &head, k, nwp);
    int age = ElapsedTime(tail,head,nwp);
    if(from<=age && (to==0 || age<to)) CHANGE_STAT[0]++; 
  }
}

/*****************
 void d_edge_ages

Sum of ages of all extant ties. This quantity changes when the clock
advances, so it needs a t_??? function.

*****************/

D_CHANGESTAT_FN(d_edge_ages_mon){
  int edgeflag, i;
  Vertex tail, head;
  
  ZERO_ALL_CHANGESTATS(i); 
  FOR_EACH_TOGGLE(i) {
    int age = ElapsedTime(tail=TAIL(i),head=HEAD(i),nwp);
    edgeflag = IS_OUTEDGE(tail, head);
    CHANGE_STAT[0] += edgeflag ? - age : 0;
    TOGGLE_IF_MORE_TO_COME(i); 
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

T_CHANGESTAT_FN(t_edge_ages_mon){
  CHANGE_STAT[0] = +N_EDGES; // Each extant edge's age increase by 1, so their sum increases by their number.
}

S_CHANGESTAT_FN(s_edge_ages_mon){
  CHANGE_STAT[0] = 0;
  for (Edge k=1; k <= N_EDGES; k++) {
    Vertex tail, head;
    FindithEdge(&tail, &head, k, nwp);
    int age = ElapsedTime(tail,head,nwp);
    CHANGE_STAT[0] += age;
  }
}


/*****************
 void d_edgecov_ages

Weighted sum of ages of all extant ties. This quantity changes when the clock
advances, so it needs a t_??? function.

*****************/

D_CHANGESTAT_FN(d_edgecov_ages_mon){
  int noffset = BIPARTITE, nrow;
  if(noffset > 0){
    nrow = noffset;
  }else{
    nrow = INPUT_PARAM[0];
  }

  int i;
  Vertex tail, head;
  
  ZERO_ALL_CHANGESTATS(i); 
  FOR_EACH_TOGGLE(i) {
    int age = ElapsedTime(tail=TAIL(i),head=HEAD(i),nwp);
    int edgeflag = IS_OUTEDGE(tail, head);
    double val = INPUT_ATTRIB[(head-1-noffset)*nrow+(tail-1)];
    CHANGE_STAT[0] += edgeflag ? - age*val : 0;
    TOGGLE_IF_MORE_TO_COME(i); 
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

T_CHANGESTAT_FN(t_edgecov_ages_mon){
  int noffset = BIPARTITE, nrow;
  if(noffset > 0){
    nrow = noffset;
  }else{
    nrow = INPUT_PARAM[0];
  }

  CHANGE_STAT[0] = 0;
  for (Edge k=1; k <= N_EDGES; k++) {
    Vertex tail, head;
    FindithEdge(&tail, &head, k, nwp);
    double val = INPUT_ATTRIB[(head-1-noffset)*nrow+(tail-1)];
    CHANGE_STAT[0] += val;
  }
  
}

S_CHANGESTAT_FN(s_edgecov_ages_mon){
  int noffset = BIPARTITE, nrow;
  if(noffset > 0){
    nrow = noffset;
  }else{
    nrow = INPUT_PARAM[0];
  }

  CHANGE_STAT[0] = 0;
  for (Edge k=1; k <= N_EDGES; k++) {
    Vertex tail, head;
    FindithEdge(&tail, &head, k, nwp);
    double val = INPUT_ATTRIB[(head-1-noffset)*nrow+(tail-1)];
    int age = ElapsedTime(tail,head,nwp);
    CHANGE_STAT[0] += age*val;
  }
}


