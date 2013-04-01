#include "changestats_homoproportion.h"

D_CHANGESTAT_FN(d_homoproportion) { 
  double change;
  int i;
  double multfactor;
  int edgeflag, num01ties, num11ties;
  Vertex ninputs, tail, head;
  Edge e;
  
  ninputs = N_INPUT_PARAMS - N_NODES - 1;
  multfactor = INPUT_PARAM[N_INPUT_PARAMS-1];
  num01ties = 0;
  num11ties = 0;
  ZERO_ALL_CHANGESTATS(i);
  /* match on attributes */
  for (tail=1; tail<=N_NODES; tail++){
   STEP_THROUGH_OUTEDGES(tail, e, head) { /* step through outedges of tail */
    if (INPUT_PARAM[tail+ninputs-1] == INPUT_PARAM[head+ninputs-1]) {
     num11ties++;
    }else{
     num01ties++;
    }
   }
  }
  
  change=0.0;
  FOR_EACH_TOGGLE(i) {
    tail = TAIL(i);
    head = HEAD(i);
    if(num11ties+num01ties>0){
     change -= num11ties*multfactor/((double)(num11ties+num01ties));
    }
    edgeflag = IS_OUTEDGE(tail, head);
    if (INPUT_PARAM[tail+ninputs-1] == INPUT_PARAM[head+ninputs-1]) {
     if(edgeflag){num11ties--;}else{num11ties++;}
    }else{
     if(edgeflag){num01ties--;}else{num01ties++;}
    }
    if(num11ties+num01ties>0){
     change += num11ties*multfactor/((double)(num11ties+num01ties));
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  CHANGE_STAT[0] = change;
  UNDO_PREVIOUS_TOGGLES(i);
}
