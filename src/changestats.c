/*  File src/changestats.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2013 Statnet Commons
 */
#include "changestats.h"

/********************  changestats:  A    ***********/
/*****************                       
 changestat: d_absdiff
*****************/
D_CHANGESTAT_FN(d_absdiff) { 
  double change, p;
  Vertex tail, head;
  int i;

  /* *** don't forget tail -> head */
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    tail = TAIL(i); 
    head = HEAD(i);
    p = INPUT_ATTRIB[0];
    if(p==1.0){
      change = fabs(INPUT_ATTRIB[tail] - INPUT_ATTRIB[head]);
    } else {
      change = pow(fabs(INPUT_ATTRIB[tail] - INPUT_ATTRIB[head]), p);
    }
    CHANGE_STAT[0] += IS_OUTEDGE(tail,head) ? -change : change;
    TOGGLE_IF_MORE_TO_COME(i); /* Needed in case of multiple toggles */
  }
  UNDO_PREVIOUS_TOGGLES(i); /* Needed on exit in case of multiple toggles */
}

/*****************
 changestat: d_absdiffcat
*****************/
D_CHANGESTAT_FN(d_absdiffcat) { 
  double change, absdiff, NAsubstitute, tailval, headval;
  Vertex tail, head, ninputs;
  int i, j;
  
  ninputs = N_INPUT_PARAMS - N_NODES;
  NAsubstitute = INPUT_PARAM[ninputs-1];
  ZERO_ALL_CHANGESTATS(i);

  /* *** don't forget tail -> head */
  FOR_EACH_TOGGLE(i) {
    change = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1.0 : 1.0;
    tailval = INPUT_ATTRIB[tail-1];
    headval = INPUT_ATTRIB[head-1];
    if (tailval == NAsubstitute ||  headval == NAsubstitute) absdiff = NAsubstitute;
    else absdiff = fabs(tailval - headval);
	  if (absdiff>0) {
      for (j=0; j<N_CHANGE_STATS; j++) {
        CHANGE_STAT[j] += (absdiff==INPUT_PARAM[j]) ? change : 0.0;
      }
    }
    TOGGLE_IF_MORE_TO_COME(i); /* Needed in case of multiple toggles */
  }
  UNDO_PREVIOUS_TOGGLES(i); /* Needed on exit in case of multiple toggles */
}

/*****************
 changestat: d_altkstar
*****************/
D_CHANGESTAT_FN(d_altkstar) { 
  int i, isedge;
  double lambda, oneexpl, change;
  Vertex tail, head, taild, headd=0;
  
  change = 0.0;
  lambda = INPUT_PARAM[0];
  oneexpl = 1.0-1.0/lambda;

  /* *** don't forget tail -> head */
  FOR_EACH_TOGGLE(i) {
    isedge = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i));
    taild = OUT_DEG[tail] + IN_DEG[tail] - isedge;
    headd = OUT_DEG[head] + IN_DEG[head] - isedge;
    if(taild!=0){
      change += (1-2*isedge)*(1.0-pow(oneexpl,(double)taild));
    }
    if(headd!=0){
      change += (1-2*isedge)*(1.0-pow(oneexpl,(double)headd));
    }
    TOGGLE_IF_MORE_TO_COME(i); /* Needed in case of multiple toggles */
  }
  CHANGE_STAT[0] = change*lambda;  
  UNDO_PREVIOUS_TOGGLES(i); /* Needed on exit in case of multiple toggles */
}

/*****************
 changestat: d_asymmetric
*****************/
D_CHANGESTAT_FN(d_asymmetric) { 
  double matchval, change;
  Vertex tail, head;
  int i, j, ninputs, noattr;

  ninputs = N_INPUT_PARAMS - N_NODES;
  noattr = (N_INPUT_PARAMS == 0);
  ZERO_ALL_CHANGESTATS(i);

  /* *** don't forget tail -> head */
  FOR_EACH_TOGGLE(i) {
    tail = TAIL(i);
    head = HEAD(i);
    change = (IS_OUTEDGE(tail, head)==IS_OUTEDGE(head, tail) ? 1.0 : -1.0) ;
    if (noattr) { /* "plain vanilla" asymmetric, without node attributes */
      CHANGE_STAT[0] += change;
    } else { /* Only consider asymmetrics where node attributes match */
      matchval = INPUT_PARAM[tail+ninputs-1];
      if (matchval == INPUT_PARAM[head+ninputs-1]) { /* We have a match! */
        if (ninputs==0) {/* diff=F in network statistic specification */
          CHANGE_STAT[0] += change;
        } else { /* diff=T */
          for (j=0; j<ninputs; j++) {
            if (matchval == INPUT_PARAM[j]) 
              CHANGE_STAT[j] += change;
          }
        }
      }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/********************  changestats:  B    ***********/
  /* For all bipartite networks:
   It is assumed that in this bipartite network, the only edges are
  of the form (b1, b2), where b1 is always strictly less
  than b2.  In other words, the degree of a b1 is equivalent
  to its outdegree and the degree of a b2 is equivalent to its
  indegree. */

/*****************
 changestat: d_b1concurrent
*****************/
D_CHANGESTAT_FN(d_b1concurrent) { 
  int i, echange;
  Vertex b1, b1deg;

  /* *** don't forget tail -> head */  
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i) {
    b1 = TAIL(i);
    echange = IS_OUTEDGE(b1, HEAD(i)) ? -1 : 1;
    b1deg = OUT_DEG[b1];
    CHANGE_STAT[0] += (b1deg + echange > 1) - (b1deg > 1);
    TOGGLE_IF_MORE_TO_COME(i); /* Needed in case of multiple toggles */
  }
  UNDO_PREVIOUS_TOGGLES(i); /* Needed on exit in case of multiple toggles */
}

/*****************
 changestat: d_b1concurrent_by_attr
*****************/
D_CHANGESTAT_FN(d_b1concurrent_by_attr) { 
  /*The inputparams are assumed to be set up as follows:
    The first 2*nstats values are in pairs:  (degree, attrvalue)
    The values following the first 2*nstats values are the nodal attributes. */
  int i, j, echange, b1attr;
  Vertex b1, b1deg;

  /* *** don't forget tail -> head */
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    b1 = TAIL(i);
    echange = IS_OUTEDGE(b1, HEAD(i)) ? -1 : 1;
    b1deg = OUT_DEG[b1];
    b1attr = INPUT_PARAM[N_CHANGE_STATS + b1 - 1]; 
    for(j = 0; j < N_CHANGE_STATS; j++) {
      if (b1attr == INPUT_PARAM[j]) { /* we have attr match */
        CHANGE_STAT[j] += (b1deg + echange > 1) - (b1deg > 1);
      }
    }
    TOGGLE_IF_MORE_TO_COME(i); /* Needed in case of multiple toggles */
  }
  UNDO_PREVIOUS_TOGGLES(i); /* Needed on exit in case of multiple toggles */
}

// A macro indicating whether x is in [from,to)
#define FROM_TO(x, from, to) ((x)>=(from) && (x)<(to))

/*****************
 changestat: d_b1degrange
*****************/
D_CHANGESTAT_FN(d_b1degrange) { 
  int i, j, echange;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    Vertex b1;
    echange=(EdgetreeSearch(b1=TAIL(i), HEAD(i), nwp->outedges)==0)? 1:-1;
    Vertex b1deg = OUT_DEG[b1];
    for(j = 0; j < N_CHANGE_STATS; j++) {
      Vertex from = INPUT_PARAM[2*j], to = INPUT_PARAM[2*j+1];
      CHANGE_STAT[j] += FROM_TO(b1deg + echange, from, to) - FROM_TO(b1deg, from, to);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}
 
/*****************
 changestat: d_b1degrange_by_attr
*****************/
D_CHANGESTAT_FN(d_b1degrange_by_attr) { 
  /* The inputparams are assumed to be set up as follows:
  The first 3*nstats values are in triples:  (from, to, attrvalue)
  The values following the first 2*nstats values are the nodal attributes.
  */
  int i, j;
  Vertex *od;
  
  od=OUT_DEG;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    Vertex b1;
    int echange = IS_OUTEDGE(b1=TAIL(i), HEAD(i)) ? -1:1;
    Vertex b1deg = od[b1];
    int b1attr = INPUT_PARAM[3*N_CHANGE_STATS + b1 - 1]; 
    for(j = 0; j < N_CHANGE_STATS; j++){
      Vertex from = INPUT_PARAM[3*j], to = INPUT_PARAM[3*j + 1];
      int testattr = INPUT_PARAM[3*j + 2]; 
      if (b1attr == testattr)  /* we have tail attr match */
        CHANGE_STAT[j] += FROM_TO(b1deg + echange, from, to) - FROM_TO(b1deg, from, to);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_b1degrange_w_homophily
*****************/
D_CHANGESTAT_FN(d_b1degrange_w_homophily) { 
  /*  The inputparams are assumed to be set up as follows:
  The first 2*nstats values are the values of b1degrange
  The values following the first 2*nstats values are the nodal attributes.
  */
  int i, j;
  double *nodeattr;
  Edge e;

  nodeattr = mtp->inputparams + N_CHANGE_STATS*2 - 1;  

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    Vertex b1=TAIL(i), b2=HEAD(i);
    int b1attr = nodeattr[b1], b2attr = nodeattr[b2];
    if (b1attr == b2attr) { /* They match; otherwise don't bother */
      int echange = IS_OUTEDGE(b1, b2) ? -1:1;
      Vertex b1deg=0, v;
      STEP_THROUGH_OUTEDGES(b1, e, v) { b1deg += (nodeattr[v]==b1attr); }
      for(j = 0; j < N_CHANGE_STATS; j++) {
        Vertex from = INPUT_PARAM[2*j], to = INPUT_PARAM[2*j+1];
        CHANGE_STAT[j] += FROM_TO(b1deg + echange, from, to) - FROM_TO(b1deg, from, to);
      }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}                                        

#undef FROM_TO

/*****************
 changestat: d_b1degree
*****************/
D_CHANGESTAT_FN(d_b1degree) { 
  int i, j, echange;
  Vertex b1, b1deg, d;

  /* *** don't forget tail -> head */  
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    b1 = TAIL(i);
    echange = IS_OUTEDGE(b1, HEAD(i)) ? -1 : 1;
    b1deg = OUT_DEG[b1];
    for(j = 0; j < N_CHANGE_STATS; j++) {
      d = (Vertex)(INPUT_PARAM[j]);
      CHANGE_STAT[j] += (b1deg + echange == d) - (b1deg == d);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_b1degree_by_attr
*****************/
D_CHANGESTAT_FN(d_b1degree_by_attr) { 
  /* The inputparams are assumed to be set up as follows:
     The first 2*nstats values are in pairs:  (degree, attrvalue)
     The values following the first 2*nstats values are the nodal attributes. */
  int i, j, echange, b1attr;
  Vertex b1, b1deg, d;
  
  /* *** don't forget tail -> head */  
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    b1 = TAIL(i);
    echange = IS_OUTEDGE(b1, HEAD(i)) ? -1 : 1;
    b1deg = OUT_DEG[b1];
    b1attr = INPUT_PARAM[2*N_CHANGE_STATS + b1 - 1]; 
    for(j = 0; j < N_CHANGE_STATS; j++) {
      if (b1attr == INPUT_PARAM[2*j+1]) { /* we have attr match */
        d = (Vertex)INPUT_PARAM[2*j];
        CHANGE_STAT[j] += (b1deg + echange == d) - (b1deg == d);
      }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_b1factor
*****************/
D_CHANGESTAT_FN(d_b1factor) { 
  double s, factorval;
  Vertex b1;
  int i, j;
  
  /* *** don't forget tail -> head */  
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    b1 = TAIL(i);
    s = IS_OUTEDGE(b1, HEAD(i)) ? -1.0 : 1.0;
    for (j=0; j<(N_CHANGE_STATS); j++) {
      factorval = (INPUT_PARAM[j]);
      CHANGE_STAT[j] += ((INPUT_ATTRIB[b1-1] != factorval) ? 0.0 : s);
    }
    TOGGLE_IF_MORE_TO_COME(i); /* Needed in case of multiple toggles */
  }
  UNDO_PREVIOUS_TOGGLES(i); /* Needed on exit in case of multiple toggles */
}

/*****************
 changestat: d_b1starmix
*****************/
D_CHANGESTAT_FN(d_b1starmix) { 
  double change;
  int edgeflag, i, j, kmo;
  Edge e;
  Vertex tail, head, node3, nnodes, taild;
  int nstats;
  double tailattr, headattr;
  
  nstats  = (int)N_CHANGE_STATS;
  nnodes = N_NODES;
  kmo = (int)INPUT_PARAM[0] - 1;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    /* edgeflag is 1 if edge exists and will disappear
    edgeflag is 0 if edge DNE and will appear */
    edgeflag = IS_OUTEDGE(tail = TAIL(i), head = HEAD(i));
    tailattr = INPUT_ATTRIB[tail-1];
    headattr = INPUT_ATTRIB[head-1];
    taild = - edgeflag; /* if edge exists set to -1 because it will be recounted */

    STEP_THROUGH_OUTEDGES(tail, e, node3) { /* step through outedges of tail */
      if(headattr == INPUT_ATTRIB[node3-1]){++taild;}
    }
    for(j=0; j < N_CHANGE_STATS; j++) {
      if (INPUT_ATTRIB[nnodes+j] == tailattr && 
      INPUT_ATTRIB[nnodes+nstats+j] == headattr) {
        change = CHOOSE(taild, kmo); 
        CHANGE_STAT[j] += (edgeflag ? - change : change); 
      }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_b1starmixhomophily
*****************/
D_CHANGESTAT_FN(d_b1starmixhomophily) { 
  double change;
  int edgeflag, i, j, kmo;
  Edge e;
  Vertex tail, head, node3, nnodes, taild;
  double tailattr, headattr;
  
  nnodes = N_NODES;
  kmo = (int)INPUT_PARAM[0] - 1;
  

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    /* edgeflag is 1 if edge exists and will disappear
    edgeflag is 0 if edge DNE and will appear */
    edgeflag = IS_OUTEDGE(tail = TAIL(i), head = HEAD(i));
    tailattr = INPUT_ATTRIB[tail-1];
    headattr = INPUT_ATTRIB[head-1];
    taild = - edgeflag; /* if edge exists set to -1 because it will be recounted */

    STEP_THROUGH_OUTEDGES(tail, e, node3) { /* step through outedges of tail */
      if(headattr == INPUT_ATTRIB[node3-1]){++taild;}
    }
    for(j=0; j < N_CHANGE_STATS; j++) {
      if (INPUT_ATTRIB[nnodes+j] == tailattr) {
        change = CHOOSE(taild, kmo); 
        CHANGE_STAT[j] += (edgeflag ? - change : change); 
      }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_b1twostar
*****************/
D_CHANGESTAT_FN(d_b1twostar) { 
  double change;
  int i, j;
  Edge e;
  Vertex tail, head, node3, nnodes;
  int nstats;
  double tailattr, headattr, n3attr;
  
  nstats  = (int)N_CHANGE_STATS;
  nnodes = N_NODES;
  
  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    change = IS_OUTEDGE(tail = TAIL(i), head = HEAD(i))? -1.0 : 1.0 ;
    tailattr = INPUT_PARAM[tail-1];
    headattr = INPUT_PARAM[head-1];

    STEP_THROUGH_OUTEDGES(tail, e, node3) { /* step through outedges of tail */
      n3attr = INPUT_PARAM[node3-1];
      for(j=0; j < N_CHANGE_STATS; j++) {
        if (node3 != head && INPUT_PARAM[nnodes + j] == tailattr && 
            INPUT_PARAM[nnodes + nstats + j] == MIN(headattr, n3attr) &&
            INPUT_PARAM[nnodes + 2*nstats + j] == MAX(headattr, n3attr)) {
          CHANGE_STAT[j] += change;
        }
      }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_b2concurrent
*****************/
D_CHANGESTAT_FN(d_b2concurrent) { 
  int i, echange;
  Vertex b2, b2deg;

  /* *** don't forget tail -> head */    
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i) {
    b2 = HEAD(i);
    echange = IS_OUTEDGE(TAIL(i), b2) ? -1 : 1;
    b2deg = IN_DEG[b2];
    CHANGE_STAT[0] += (b2deg + echange > 1) - (b2deg > 1);
    TOGGLE_IF_MORE_TO_COME(i); /* Needed in case of multiple toggles */
  }
  UNDO_PREVIOUS_TOGGLES(i); /* Needed on exit in case of multiple toggles */
}

/*****************
 changestat: d_b2concurrent_by_attr
*****************/
D_CHANGESTAT_FN(d_b2concurrent_by_attr) { 
  /*The inputparams are assumed to be set up as follows:
    The first 2*nstats values are in pairs:  (degree, attrvalue)
    The values following the first 2*nstats values are the nodal attributes.*/
  int i, j, echange, b2attr;
  Vertex b2, b2deg;
  

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    b2 = HEAD(i);
    echange = IS_OUTEDGE(TAIL(i), b2) ? -1 : 1;
    b2deg = IN_DEG[b2];
    b2attr = INPUT_PARAM[N_CHANGE_STATS + b2 - 1 - BIPARTITE];
    for(j = 0; j < N_CHANGE_STATS; j++) {
      if (b2attr == INPUT_PARAM[j]) { /* we have attr match */
        CHANGE_STAT[j] += (b2deg + echange > 1) - (b2deg > 1);
      }
    }
    TOGGLE_IF_MORE_TO_COME(i); /* Needed in case of multiple toggles */
  }
  UNDO_PREVIOUS_TOGGLES(i); /* Needed on exit in case of multiple toggles */
}

// A macro indicating whether x is in [from,to)
#define FROM_TO(x, from, to) ((x)>=(from) && (x)<(to))

/*****************
 changestat: d_b2degrange
*****************/
D_CHANGESTAT_FN(d_b2degrange) { 
  int i, j, echange;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    Vertex b2;
    echange=(EdgetreeSearch(TAIL(i), b2=HEAD(i), nwp->outedges)==0)? 1:-1;
    Vertex b2deg = IN_DEG[b2];
    for(j = 0; j < N_CHANGE_STATS; j++) {
      Vertex from = INPUT_PARAM[2*j], to = INPUT_PARAM[2*j+1];
      CHANGE_STAT[j] += FROM_TO(b2deg + echange, from, to) - FROM_TO(b2deg, from, to);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}
 
/*****************
 changestat: d_b2degrange_by_attr
*****************/
D_CHANGESTAT_FN(d_b2degrange_by_attr) { 
  /* The inputparams are assumed to be set up as follows:
  The first 3*nstats values are in triples:  (from, to, attrvalue)
  The values following the first 2*nstats values are the nodal attributes.
  */
  int i, j;
  Vertex *id;
  
  id=IN_DEG;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    Vertex b2;
    int echange = IS_OUTEDGE(TAIL(i), b2=HEAD(i)) ? -1:1;
    Vertex b2deg = id[b2];
    int b1attr = INPUT_PARAM[3*N_CHANGE_STATS + b2 - 1 - BIPARTITE]; 
    for(j = 0; j < N_CHANGE_STATS; j++){
      Vertex from = INPUT_PARAM[3*j], to = INPUT_PARAM[3*j + 1];
      int testattr = INPUT_PARAM[3*j + 2]; 
      if (b1attr == testattr)  /* we have tail attr match */
        CHANGE_STAT[j] += FROM_TO(b2deg + echange, from, to) - FROM_TO(b2deg, from, to);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_b2degrange_w_homophily
*****************/
D_CHANGESTAT_FN(d_b2degrange_w_homophily) { 
  /*  The inputparams are assumed to be set up as follows:
  The first 2*nstats values are the values of b2degrange
  The values following the first 2*nstats values are the nodal attributes.
  */
  int i, j;
  double *nodeattr;
  Edge e;

  nodeattr = mtp->inputparams + N_CHANGE_STATS*2 - 1;  

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    Vertex b1=TAIL(i), b2=HEAD(i);
    int b1attr = nodeattr[b1], b2attr = nodeattr[b2];
    if (b1attr == b2attr) { /* They match; otherwise don't bother */
      int echange = IS_OUTEDGE(b1, b2) ? -1:1;
      Vertex b2deg=0, v;
      STEP_THROUGH_INEDGES(b2, e, v) { b2deg += (nodeattr[v]==b1attr); }
      for(j = 0; j < N_CHANGE_STATS; j++) {
        Vertex from = INPUT_PARAM[2*j], to = INPUT_PARAM[2*j+1];
        CHANGE_STAT[j] += FROM_TO(b2deg + echange, from, to) - FROM_TO(b2deg, from, to);
      }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}                                        

#undef FROM_TO


/*****************
 changestat: d_b2degree
*****************/
D_CHANGESTAT_FN(d_b2degree) { 
  int i, j, echange;
  Vertex b2, b2deg, d;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    b2 = HEAD(i);
    echange = IS_OUTEDGE(TAIL(i), b2) ? -1 : 1;
    b2deg = IN_DEG[b2];
    for(j = 0; j < N_CHANGE_STATS; j++) {
      d = (Vertex)(INPUT_PARAM[j]);
      CHANGE_STAT[j] += (b2deg + echange == d) - (b2deg == d);
    }
    TOGGLE_IF_MORE_TO_COME(i); /* Needed in case of multiple toggles */
  }
  UNDO_PREVIOUS_TOGGLES(i); /* Needed on exit in case of multiple toggles */
}

/*****************
 changestat: d_b2degree_by_attr
*****************/
D_CHANGESTAT_FN(d_b2degree_by_attr) { 
  /*The inputparams are assumed to be set up as follows:
    The first 2*nstats values are in pairs:  (degree, attrvalue)
    The values following the first 2*nstats values are the nodal attributes. */
  int i, j, echange, b2attr;
  Vertex b2, b2deg, d;
  

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    b2 = HEAD(i);
    echange = IS_OUTEDGE(TAIL(i), b2) ? -1 : 1;
    b2deg = IN_DEG[b2];
    b2attr = INPUT_PARAM[2*N_CHANGE_STATS + b2 - 1 - BIPARTITE];
    for(j = 0; j < N_CHANGE_STATS; j++) {
      if (b2attr == INPUT_PARAM[2*j+1]) { /* we have attr match */
        d = (Vertex)INPUT_PARAM[2*j];
        CHANGE_STAT[j] += (b2deg + echange == d) - (b2deg == d);
      }
    }
    TOGGLE_IF_MORE_TO_COME(i); /* Needed in case of multiple toggles */
  }
  UNDO_PREVIOUS_TOGGLES(i); /* Needed on exit in case of multiple toggles */
}

/*****************
 changestat: d_b2factor
*****************/
D_CHANGESTAT_FN(d_b2factor) { 
  double s, factorval;
  Vertex nb1, b2;
  int i, j;
  

  /* *** don't forget tail -> head */    
  nb1 = BIPARTITE;
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    b2 = HEAD(i);
    s = IS_OUTEDGE(TAIL(i), b2) ? -1.0 : 1.0;
    for (j=0; j<(N_CHANGE_STATS); j++) {
      factorval = (INPUT_PARAM[j]);
      CHANGE_STAT[j] += ((INPUT_ATTRIB[b2-nb1-1] != factorval) ? 0.0 : s);
    }
    TOGGLE_IF_MORE_TO_COME(i); /* Needed in case of multiple toggles */
  }
  UNDO_PREVIOUS_TOGGLES(i); /* Needed on exit in case of multiple toggles */
}

/*****************
 changestat: d_b2starmix
*****************/
D_CHANGESTAT_FN(d_b2starmix) { 
  double change;
  int edgeflag, i, j, kmo;
  Edge e;
  Vertex tail, head, node3, nnodes, headd;
  int nstats;
  double tailattr, headattr;
  
  nstats  = (int)N_CHANGE_STATS;
  nnodes = N_NODES;
  kmo = (int)INPUT_PARAM[0] - 1;
  

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    /* edgeflag is 1 if edge exists and will disappear
    edgeflag is 0 if edge DNE and will appear */
    edgeflag = IS_OUTEDGE(tail = TAIL(i), head = HEAD(i));
    tailattr = INPUT_ATTRIB[tail-1];
    headattr = INPUT_ATTRIB[head-1];
    headd = - edgeflag; /* if edge exists set to -1 because it will be recounted */

    STEP_THROUGH_INEDGES(head, e, node3) { /* step through inedges of head */
      if(tailattr == INPUT_ATTRIB[node3-1]){++headd;}
    }
    for(j=0; j < N_CHANGE_STATS; j++) {
      if (INPUT_ATTRIB[nnodes+j] == tailattr && 
      INPUT_ATTRIB[nnodes+nstats+j] == headattr) {
        change = CHOOSE(headd, kmo); 
        CHANGE_STAT[j] += (edgeflag ? - change : change); 
      }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_b2starmixhomophily
*****************/
D_CHANGESTAT_FN(d_b2starmixhomophily) { 
  double change;
  int edgeflag, i, j, kmo;
  Edge e;
  Vertex tail, head, node3, nnodes, headd;
  double tailattr, headattr;
  
  nnodes = N_NODES;
  kmo = (int)INPUT_PARAM[0] - 1;
  

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    /* edgeflag is 1 if edge exists and will disappear
    edgeflag is 0 if edge DNE and will appear */
    edgeflag = IS_OUTEDGE(tail = TAIL(i), head = HEAD(i));
    tailattr = INPUT_ATTRIB[tail-1];
    headattr = INPUT_ATTRIB[head-1];
    headd = - edgeflag; /* if edge exists set to -1 because it will be recounted */

    STEP_THROUGH_INEDGES(head, e, node3) { /* step through inedges of head */
      if(tailattr == INPUT_ATTRIB[node3-1]){++headd;}
    }
    for(j=0; j < N_CHANGE_STATS; j++) {
      if (INPUT_ATTRIB[nnodes+j] == headattr) {
        change = CHOOSE(headd, kmo); 
        CHANGE_STAT[j] += (edgeflag ? - change : change); 
      }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_b2twostar
*****************/
D_CHANGESTAT_FN(d_b2twostar) { 
  double change;
  int i, j;
  Edge e;
  Vertex tail, head, node3, nnodes;
  int nstats;
  double tailattr, headattr, n3attr;
  
  nstats  = (int)N_CHANGE_STATS;
  nnodes = N_NODES;
  

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    change = IS_OUTEDGE(tail = TAIL(i), head = HEAD(i))? -1.0 : 1.0 ;
    tailattr = INPUT_PARAM[tail-1];
    headattr = INPUT_PARAM[head-1];

    STEP_THROUGH_INEDGES(head, e, node3) { /* step through inedges of head */
      n3attr = INPUT_PARAM[node3-1];
      for(j=0; j < N_CHANGE_STATS; j++) {
        if (node3 != tail && INPUT_PARAM[nnodes + j] == headattr && 
            INPUT_PARAM[nnodes + nstats + j] == MIN(tailattr, n3attr) &&
            INPUT_PARAM[nnodes + 2*nstats + j] == MAX(tailattr, n3attr)) {
          CHANGE_STAT[j] += change;
        }
      }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_balance
*****************/
D_CHANGESTAT_FN(d_balance) { 
  int i, edgeflag, a, b, c, d, e, edgecount, t300, 
  t210, t120C, t120U, t120D, t201, t030C, t030T, t111U, 
  t111D, t021C, t021U, t021D, t102, t012; /* , t003; */
  Vertex node3, tail, head;


  /* *** don't forget tail -> head */    
  CHANGE_STAT[0] = 0.0;
  if (DIRECTED) { /* directed version */
    FOR_EACH_TOGGLE(i) {
      tail = TAIL(i);
      head = HEAD(i);
      edgeflag = IS_OUTEDGE(tail, head);
      t300 = 0;
      t210 = 0;
      t120C = 0;  t120U = 0;   t120D = 0;  t201 = 0;
      t030C = 0;  t030T = 0;   t111U = 0;  t111D = 0;
      t021C = 0;  t021U = 0;   t021D = 0;  t102 = 0;
      t012 = 0;

      if (MIN_OUTEDGE(head)!=0 || MIN_INEDGE(head)!=0 || 
      MIN_OUTEDGE(tail)!=0 || MIN_INEDGE(tail)!=0) {

        /* ****** loop through node3 ****** */
          for (node3=1; node3 <= N_NODES; node3++) { 
            if (node3 != tail && node3 != head) {  
              a = IS_OUTEDGE(head, tail);
              b = IS_OUTEDGE(head, node3);
              c = IS_OUTEDGE(node3, head);
              d = IS_OUTEDGE(node3, tail);
              e = IS_OUTEDGE(tail, node3);
              edgecount = (a + b + c + d + e);
              
              switch(edgecount) {  
                case 0:   /* 012 */
                ++t012;
                
                case 1: {  /* 021C, 021U, 021D, 102 */
                  if ((b == 1) || (d == 1))
                    ++t021C;
                  if (c == 1)
                    ++t021U;
                  if (e == 1)
                    ++t021D;
                  if (a == 1)
                    ++t102;
                }
                break;

                case 2:  { /* 030C, 030T, 111U, 111D */
                  if ((b + d) == 2)       
                    ++t030C;
                  if (((b + e) == 2) || ((c + d) == 2) || ((c + e) == 2))
                    ++t030T;
                  if (((a + b) == 2) || ((a + e) == 2) || ((d + e) == 2))
                    ++t111U;
                  if (((a + c) == 2) || ((a + d) == 2) || ((b + c) == 2))
                    ++t111D;
                }
                break;
            
                case 3: {   /* 120C, 120U, 120D, 201 */
                  if (a == 1) {
                    if (((b + d) == 2) || ((c + e) == 2))
                      ++t120C;
                    if ((b + e) == 2)
                      ++t120U;
                    if ((c + d) == 2)
                      ++t120D;
                    if (((b + c) == 2) || ((d + e) == 2))
                      ++t201;
                  } else { 
                    if (b == 1) {
                      if (((c + d) == 2) || ((d + e) == 2))
                        ++t120C;
                      if ((c + e) == 2)
                        ++t120D;
                    } else {
                      ++t120U;
                    }
                  } 
                }
                break;
             
                case 4:   /* 210 */
                ++t210;
                break;
             
                case 5:   /* 300 */            
                ++t300;
                break;
              }

              switch(edgecount) {  
                case 1:   /* 102, 021D, 021U, 021C */
                --t012;
                break;
                
                case 2: {   /* 030C, 030T, 111U, 111D */ 
                  if (((a + c) == 2) || ((a + e) == 2) || ((b + d) == 2) || 
                    ((c + e) == 2)) 
                  --t021C;
                  if (((a + d) == 2) || ((b + e) == 2))
                    --t021U;
                  if (((a + b) == 2) || ((c + d) == 2))
                    --t021D;
                  if (((b + c) == 2) || ((d + e) == 2))
                    --t102;
                } 
                break;

                case 3: {  /* 201, 120D, 120U, 120C */
                  if (a == 1) {
                    if ((c + e) == 2)       
                      --t030C;
                    if (((c + d) == 2) || ((b + e) == 2) || ((b + d) == 2))
                      --t030T;
                    if ((b + c) == 2)
                      --t111U;
                    if ((d + e) == 2)
                      --t111D;
                  } else {
                    if (b == 1) {
                      if ((c + d) == 2)
                        --t111U;
                      if (((c + e) == 2) || ((d + e) == 2))
                        --t111D;
                    }
                    else
                      --t111U;
                  }
                }
                break;
          
                case 4: {   /* 210 */
                  if (a == 1) {
                    if (((b + c + e) == 3) || ((c + d + e) == 3))
                      --t120C;
                    if ((b + c + d) == 3)
                      --t120U;
                    if ((b + d + e) == 3)
                      --t120D;
                  } else {
                    if ((b + c + d + e) == 4)
                      --t201;
                  } 
                }
                break;
                
                case 5:   /* 300 */            
                --t210;
                break;
              }
            }
          }    /* ******  move to next node3 ******** */
      }
      else 
        t012 = t012 + (N_NODES - 2);  
      
      /*        t003 = (t300+t210+t120C+t120U+t120D+t201+t030C+t030T); 
      t003 = t003+(t111U+t111D+t021C+t021U+t021D+t102+t012); */
      b = t102 + t300; 
      CHANGE_STAT[0] += edgeflag ? -(double)b : (double)b;

      TOGGLE_IF_MORE_TO_COME(i);
    }

  /* *** don't forget tail -> head */    
  }else{ /*  undirected */
    FOR_EACH_TOGGLE(i) {
      tail = TAIL(i); 
      head = HEAD(i);
      edgeflag = IS_OUTEDGE(tail, head);
      t300 = 0; t201 = 0; t102 = 0; t012 = 0;
      
      if (MIN_OUTEDGE(head)!=0 || MIN_INEDGE(head)!=0 ||
      MIN_OUTEDGE(tail)!=0 || MIN_INEDGE(tail)!=0) {
        
        /* ****** loop through node3 ****** */
        for (node3=1; node3 <= N_NODES; node3++) { 
          if (node3 != tail && node3 != head) {
            a = IS_UNDIRECTED_EDGE(node3, head);
            b = IS_UNDIRECTED_EDGE(node3, tail);
            edgecount = (a + b);
            
            switch(edgecount){  
              case 0: {  /* 012 */
                ++t102;
                --t012;
              }
              break;
              
              case 1: {  /* 021C, 021U, 021D, 102 */
                ++t201;
                --t102;
              }
              break;
              
              case 2: {  /* 030C, 030T, 111U, 111D */
                ++t300;
                --t201;
              }
              break;
            }
          }
          
        }    /* ******  move to next node3 ******** */
      } else 
      t102 = t102 + (N_NODES - 2);  
      
      /* t003 = (t102+t201+t300); */
      b = t102 + t300; 
      CHANGE_STAT[0] += edgeflag ? -(double)b : (double)b;
      
      TOGGLE_IF_MORE_TO_COME(i);
    } /* i loop */
  }
  UNDO_PREVIOUS_TOGGLES(i); /* Needed on exit in case of multiple toggles */
}
  
/*****************
 changestat: d_boundeddegree
*****************/
D_CHANGESTAT_FN(d_boundeddegree) { 
  int i, j, echange;
  Vertex tail, head, taild, headd=0, deg;
  int nstats = (int)N_CHANGE_STATS;
  Vertex bound = (Vertex)INPUT_PARAM[nstats-1];
  
  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    tail = TAIL(i);
    head = HEAD(i);
    echange = IS_OUTEDGE(tail, head) ? -1 : 1;
    taild = OUT_DEG[tail] + IN_DEG[tail];
    headd = OUT_DEG[head] + IN_DEG[head];
    for(j = 0; j+1 < nstats; j++)	{
      deg = (Vertex)INPUT_PARAM[j];
      CHANGE_STAT[j] += (taild + echange == deg) - (taild == deg);
      CHANGE_STAT[j] += (headd + echange == deg) - (headd == deg);
    }
    CHANGE_STAT[nstats-1] += (taild + echange >= bound) - (taild >= bound);
    CHANGE_STAT[nstats-1] += (headd + echange >= bound) - (headd >= bound);    
    TOGGLE_IF_MORE_TO_COME(i); /* Needed in case of multiple toggles */
  }
  UNDO_PREVIOUS_TOGGLES(i); /* Needed on exit in case of multiple toggles */
}

/*****************
 changestat: d_boundedidegree
*****************/
D_CHANGESTAT_FN(d_boundedidegree) { 
  int i, j, echange;
  Vertex tail, taild=0, deg;
  int nstats = (int)N_CHANGE_STATS;
  Vertex bound = (Vertex)INPUT_PARAM[nstats-1];
  
  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    tail = TAIL(i);
    echange = IS_OUTEDGE(tail, HEAD(i)) ? -1 : 1;
    taild = IN_DEG[tail];
    for(j = 0; j < N_CHANGE_STATS; j++) {
      deg = (Vertex)INPUT_PARAM[j];
      CHANGE_STAT[j] += (taild + echange == deg) - (taild == deg);
    }
    CHANGE_STAT[nstats-1] += (taild + echange >= bound) - (taild >= bound);
    TOGGLE_IF_MORE_TO_COME(i); /* Needed in case of multiple toggles */
  }
  UNDO_PREVIOUS_TOGGLES(i); /* Needed on exit in case of multiple toggles */
}

/*****************
 changestat: d_boundedistar
*****************/
D_CHANGESTAT_FN(d_boundedistar) { 
  double change, headod;
  double newheadod;
  int edgeflag, i, j, k, bound;
  int p = N_CHANGE_STATS;
  Vertex head;
  
  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    /* is there an edge for this toggle */
    head = HEAD(i);
    edgeflag = IS_OUTEDGE(TAIL(i), head);
    headod = IN_DEG[head];
    newheadod = headod + (edgeflag ? -1 : 1);
    for(j=0; j < p; j++) {
      k =  ((int)INPUT_PARAM[j]);
      bound = (int)INPUT_PARAM[j+p];
      change = MIN(bound,CHOOSE(newheadod, k))-MIN(bound,CHOOSE(headod, k));
      CHANGE_STAT[j] += change;
    }
    TOGGLE_IF_MORE_TO_COME(i); /* Needed in case of multiple toggles */
  }
  UNDO_PREVIOUS_TOGGLES(i); /* Needed on exit in case of multiple toggles */
}

/*****************
 changestat: d_boundedkstar
*****************/
D_CHANGESTAT_FN(d_boundedkstar) { 
  double change, tailod, headod;
  double newtailod, newheadod;
  int edgeflag, i, j, k, bound;
  int p = N_CHANGE_STATS;
  Vertex tail, head;
  
  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    /* is there an edge for this toggle */
    tail = TAIL(i);
    head = HEAD(i);
    edgeflag = IS_OUTEDGE(tail, head);
    tailod = OUT_DEG[tail] + IN_DEG[tail];
    newtailod = tailod + (edgeflag ? -1 : 1);
    headod = OUT_DEG[head] + IN_DEG[head];
    newheadod = headod + (edgeflag ? -1 : 1);
    for(j=0; j < p; j++) {
      k =  ((int)INPUT_PARAM[j]);
      bound = (int)INPUT_PARAM[j+p];
      change = (MIN(bound,CHOOSE(newtailod, k))-MIN(bound,CHOOSE(tailod, k))) +
      (MIN(bound,CHOOSE(newheadod, k))-MIN(bound,CHOOSE(headod, k)));
      
      CHANGE_STAT[j] += change; /* (edgeflag ? - change : change); */
    }
    TOGGLE_IF_MORE_TO_COME(i); /* Needed in case of multiple toggles */
  }
  UNDO_PREVIOUS_TOGGLES(i); /* Needed on exit in case of multiple toggles */
}

/*****************
 changestat: d_boundedodegree
*****************/
D_CHANGESTAT_FN(d_boundedodegree) { 
  int i, j, echange;
  Vertex tail, taild=0, deg;
  int nstats = (int)N_CHANGE_STATS;
  Vertex bound = (Vertex)INPUT_PARAM[nstats-1];
  
  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    tail = TAIL(i);
    echange = IS_OUTEDGE(tail, HEAD(i)) ? -1 : 1;
    taild = OUT_DEG[tail];
    for(j = 0; j < N_CHANGE_STATS; j++)  {
      deg = (Vertex)INPUT_PARAM[j];
      CHANGE_STAT[j] += (taild + echange == deg) - (taild == deg);
    }
    CHANGE_STAT[nstats-1] += (taild + echange >= bound) - (taild >= bound);
    TOGGLE_IF_MORE_TO_COME(i); /* Needed in case of multiple toggles */
  }
  UNDO_PREVIOUS_TOGGLES(i); /* Needed on exit in case of multiple toggles */
}

/*****************
 changestat: d_boundedostar
*****************/
D_CHANGESTAT_FN(d_boundedostar) { 
  double change, tailod;
  double newtailod;
  int edgeflag, i, j, k, bound;
  int p = N_CHANGE_STATS;
  Vertex tail;
  
  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    /* is there an edge for this toggle */
    tail = TAIL(i);
    edgeflag = IS_OUTEDGE(tail, HEAD(i));
    tailod = OUT_DEG[tail];
    newtailod = tailod + (edgeflag ? -1 : 1);
      for(j=0; j < p; j++) {
        k =  ((int)INPUT_PARAM[j]);
        bound = (int)INPUT_PARAM[j+p];
        change = MIN(bound,CHOOSE(newtailod, k))-MIN(bound,CHOOSE(tailod, k));
        CHANGE_STAT[j] += change;
      }
      TOGGLE_IF_MORE_TO_COME(i); /* Needed in case of multiple toggles */
  }
  UNDO_PREVIOUS_TOGGLES(i); /* Needed on exit in case of multiple toggles */
}

/*****************
 changestat: d_boundedtriangle
*****************/
D_CHANGESTAT_FN(d_boundedtriangle) { 
  Edge e;
  Vertex tail, head, node3;
  double boundedchange, htcount;
  Vertex tailtri, headtri;
  int edgeflag, i;
  int bound = (int)INPUT_PARAM[0];

  /* *** don't forget tail -> head */    
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i) {
    tail = TAIL(i);
    head = HEAD(i);
    edgeflag = IS_OUTEDGE(tail, head);
    tailtri=0;
    headtri=0;
    STEP_THROUGH_OUTEDGES(tail, e, node3) {
      tailtri += CountTriangles(tail, node3, 1, 1, nwp);
    }
    STEP_THROUGH_INEDGES(tail, e, node3) {
      tailtri += CountTriangles(tail, node3, 1, 1, nwp);
	  }
    STEP_THROUGH_OUTEDGES(head, e, node3) {
      headtri += CountTriangles(head, node3, 1, 1, nwp);
    }
    STEP_THROUGH_INEDGES(head, e, node3) {
      headtri += CountTriangles(head, node3, 1, 1, nwp);
	  }
    tailtri = tailtri/2;
    headtri = headtri/2;
    htcount = CountTriangles(tail, head, 1, 1, nwp);
    boundedchange = (MIN(headtri+(edgeflag ? -1:1)*htcount,bound)-MIN(headtri,bound)+
                    MIN(tailtri+(edgeflag ? -1:1)*htcount,bound)-MIN(tailtri,bound));
    CHANGE_STAT[0] += boundedchange;
    TOGGLE_IF_MORE_TO_COME(i); /* Needed in case of multiple toggles */
  }
  UNDO_PREVIOUS_TOGGLES(i); /* Needed on exit in case of multiple toggles */
}

/*****************
 CountTriangles: called by d_boundedtriangle
*****************/
Vertex CountTriangles (Vertex tail, Vertex head, int outcount, int incount, 
		       Network *nwp) {
  Edge e;
  Vertex change;
  Vertex k;
  
  /* *** don't forget tail -> head */    
  change=0;
  if(outcount){
    for(e = EdgetreeMinimum(nwp->outedges, head);
	(k = nwp->outedges[e].value) != 0;
	e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of head */
      {
	if (EdgetreeSearch(MIN(k,tail), MAX(k,tail), nwp->outedges) != 0)
	  ++change;
      }
  }
  
  if(incount){
    for(e = EdgetreeMinimum(nwp->inedges, head); 
	(k = nwp->inedges[e].value) != 0;
	e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of head */
      {
	if (EdgetreeSearch(MIN(k,tail), MAX(k,tail), nwp->outedges) != 0)
	  ++change;
      }
  }
  return(change);
}



/********************  changestats:  C    ***********/
/*****************
 changestat: d_concurrent
*****************/
D_CHANGESTAT_FN(d_concurrent) { 
  int i, echange;
  Vertex tail, head, taildeg, headdeg;

  /* *** don't forget tail -> head */    
  CHANGE_STAT[0] = 0.0;  
  FOR_EACH_TOGGLE(i) {
    tail = TAIL(i);
    head = HEAD(i);
    echange = IS_OUTEDGE(tail, head) ? -1 : 1;
    taildeg = OUT_DEG[tail];
    headdeg = IN_DEG[head];
    if(!DIRECTED){
      taildeg += IN_DEG[tail];
      headdeg += OUT_DEG[head];
    }
    CHANGE_STAT[0] += (taildeg + echange > 1) - (taildeg > 1);
    CHANGE_STAT[0] += (headdeg + echange > 1) - (headdeg > 1);
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_concurrent_by_attr
*****************/
D_CHANGESTAT_FN(d_concurrent_by_attr) { 
  /*The inputparams are assumed to be set up as follows:
    The first 2*nstats values are in pairs:  (degree, attrvalue)
    The values following the first 2*nstats values are the nodal attributes.
  */
  int i, j, echange, tailattr, headattr;
  Vertex tail, head, taildeg, headdeg;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    tail = TAIL(i);
    head = HEAD(i);
    echange = IS_OUTEDGE(tail, head) ? -1 : 1;
    taildeg = OUT_DEG[tail];
    headdeg = IN_DEG[head];
    if(!DIRECTED){
      taildeg += IN_DEG[tail];
      headdeg += OUT_DEG[head];
    }
    tailattr = INPUT_PARAM[N_CHANGE_STATS + tail - 1]; 
    headattr = INPUT_PARAM[N_CHANGE_STATS + head - 1];
    for(j = 0; j < N_CHANGE_STATS; j++) {
      if (tailattr == INPUT_PARAM[j]) { /* we have attr match */
        CHANGE_STAT[j] += (taildeg + echange > 1) - (taildeg > 1);
      }
      if (headattr == INPUT_PARAM[j]) { /* we have attr match */
        CHANGE_STAT[j] += (headdeg + echange > 1) - (headdeg > 1);
      }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_ctriple
*****************/
D_CHANGESTAT_FN(d_ctriple) { 
  Edge e;
  Vertex tail, head, change, node3;
  int i, j;
  double tailattr, edgemult;
  
  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    tail = TAIL(i);
    head = HEAD(i);
    edgemult = IS_OUTEDGE(tail, head) ? -1.0 : 1.0;
    change = 0;
    if(N_INPUT_PARAMS > 0){ /* match on attributes */
      tailattr = INPUT_ATTRIB[tail-1];
      if(tailattr == INPUT_ATTRIB[head-1]) {
        STEP_THROUGH_OUTEDGES(head, e, node3) { /* step through outedges of head */
          if(tailattr == INPUT_ATTRIB[node3-1])
            change += IS_OUTEDGE(node3, tail);
        }
        if(N_CHANGE_STATS > 1) { /* diff = TRUE; matches must be tabled */
          for (j=0; j<N_CHANGE_STATS; j++){
            if (tailattr == INPUT_PARAM[j])
              CHANGE_STAT[j] += edgemult * change;
          }
        } else { /* diff = FALSE; all matches equivalent */
          CHANGE_STAT[0] += edgemult * change;          
        }
      }
    }else{ /* no attribute matching */
      STEP_THROUGH_OUTEDGES(head, e, node3) { /* step through outedges of head */
        change += IS_OUTEDGE(node3, tail);
      }
      CHANGE_STAT[0] += edgemult * change;
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_cycle
*****************/
D_CHANGESTAT_FN(d_cycle) { 
  int i,j,k;
  Vertex tail, head;
  long int maxlen;
  double *countv,emult;
  
  /*Perform initial setup*/
  maxlen=(long int)(INPUT_PARAM[0]);
  countv=(double *)R_alloc(sizeof(double),maxlen-1);

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    for(j=0;j<maxlen-1;j++)  /*Clear out the count vector*/
      countv[j]=0.0;
    tail = TAIL(i);
    head = HEAD(i);
    /*Count the cycles associated with this edge*/
    /* OLD COMMENTS */
    /*     Note: the ergm toggle system gets heads and tails reversed!*/
    /* NEW COMMENTS */ 
    /*     *** with the h/t swap, the <edgewise_cycle_census> function
           seems correct as written */
    /*edgewise_cycle_census(g,TAIL(i),HEAD(i),countv,maxlen,directed);*/
    edgewise_cycle_census(nwp,tail,head,countv,maxlen);

    /*Make the change, as needed*/
    /* I did not swap h/t in the comment below */
    /*edgeflag = (EdgetreeSearch(tail=TAIL(i), head=HEAD(i), g.outedges) != 0);*/
    if((!DIRECTED)&&(tail>head))
      emult = IS_OUTEDGE(head, tail) ? -1.0 : 1.0;
    else
      emult = IS_OUTEDGE(tail, head) ? -1.0 : 1.0;
    k=0;
    for(j=0;j<maxlen-1;j++)
      if(INPUT_PARAM[1+j]>0.0)
        CHANGE_STAT[k++]+=emult*countv[j];
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 edgewise_path_recurse:  Called by d_cycle
*****************/
void edgewise_path_recurse(Network *nwp, Vertex dest, Vertex curnode, 
                   Vertex *availnodes, long int availcount, long int curlen, 
                   double *countv, long int maxlen) {
  Vertex *newavail,i,j;
  long int newavailcount;
  int rflag;
  
  /*If we've found a path to the destination, increment the census vector*/ 
  if(DIRECTED||(curnode<dest)) countv[curlen] += IS_OUTEDGE(curnode, dest);
  else countv[curlen] += IS_OUTEDGE(dest, curnode);
  
  /*If possible, keep searching for novel paths*/
  if((availcount>0)&&(curlen<maxlen-2)){
    if(availcount>1){    /*Remove the current node from the available list*/
      if((newavail=(Vertex *)malloc(sizeof(Vertex)*(availcount-1)))==NULL){
        Rprintf("Unable to allocate %d bytes for available node list in edgewise_path_recurse.  Trying to terminate recursion gracefully, but your path count is probably wrong.\n",sizeof(Vertex)*(availcount-1));
        return;
      }
      j=0;
      for(i=0;i<availcount;i++)      /*Create the reduced list, fur passin'*/
        if(availnodes[i]!=curnode)
          newavail[j++]=availnodes[i];
    }else
      newavail=NULL;                 /*Set to NULL if we're out of nodes*/
    newavailcount=availcount-1;      /*Decrement the available count*/

    /*Recurse on all available nodes*/
    for(i=0;i<newavailcount;i++) {
      rflag = DIRECTED || (curnode<newavail[i]) ? 
              IS_OUTEDGE(curnode,newavail[i]) : IS_OUTEDGE(newavail[i],curnode);
      if(rflag)
        edgewise_path_recurse(nwp,dest,newavail[i],newavail,newavailcount,
                              curlen+1,countv,maxlen);
    }
    /*Free the available node list*/
    if(newavail!=NULL)
      free((void *)newavail);
  }

  /*Check for interrupts (if recursion is taking way too long...)*/
  R_CheckUserInterrupt();
}

/*****************
 edgewise_cycle_census:  Called by d_cycle
*****************/
/* *** I did NOT swap heads and tails in this function, since it
   appears to have been written with tail -> head in mind */ 

void edgewise_cycle_census(Network *nwp, Vertex tail, Vertex head, 
                           double *countv, long int maxlen) {
  /* *** don't forget tail -> head */    
  long int n,i,j;
  Vertex *availnodes;
  int rflag;

  /*Set things up*/
  n=N_NODES;

  /*First, check for a 2-cycle (but only if directed)*/
  if(DIRECTED && IS_OUTEDGE(head,tail))
    countv[0]++;
  if(n==2)
    return;                 /*Failsafe for graphs of order 2*/
  
  /*Perform the recursive path count*/
  if((availnodes=(Vertex *)malloc(sizeof(Vertex)*(n-2)))==NULL){
    Rprintf("Unable to allocate %d bytes for available node list in edgewise_cycle_census.  Exiting.\n",sizeof(Vertex)*(n-2));
    return;
  }
  j=0;                             /*Initialize the list of available nodes*/
  for(i=1;i<=n;i++)
    if((i!=head)&&(i!=tail))
      availnodes[j++]=i;
  for(i=0;i<n-2;i++) {             /*Recurse on each available vertex*/
    rflag = DIRECTED || (head < availnodes[i]) ? 
            IS_OUTEDGE(head, availnodes[i]) : IS_OUTEDGE(availnodes[i], head);
    if(rflag)
      edgewise_path_recurse(nwp,tail,availnodes[i],availnodes,n-2,1,countv,maxlen);
  }
  free((void *)availnodes);  /*Free the available node list*/
}

/********************  changestats:  D    ***********/
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
 changestat: d_degcrossprod
*****************/
D_CHANGESTAT_FN(d_degcrossprod) { 
  int i, echange;
  Vertex tail, head, taildeg, headdeg, node3;
  Edge e;
  double nedges;

  nedges = INPUT_PARAM[0];

  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i) {
    tail = TAIL(i);
    head = HEAD(i);
    echange = IS_OUTEDGE(tail, head) ? -1 : 1;
    taildeg = OUT_DEG[tail] + IN_DEG[tail];
    headdeg = OUT_DEG[head] + IN_DEG[head];
    if(echange==1){
     CHANGE_STAT[0] += (taildeg + 1)*(headdeg + 1);
     STEP_THROUGH_OUTEDGES(head, e, node3) { /* step through outedges of head */
      if(node3!=tail) CHANGE_STAT[0] += (OUT_DEG[node3] + IN_DEG[node3]);
     }
     STEP_THROUGH_INEDGES(head, e, node3) { /* step through inedges of head */
      if(node3!=tail) CHANGE_STAT[0] += (OUT_DEG[node3] + IN_DEG[node3]);
     }
     STEP_THROUGH_OUTEDGES(tail, e, node3) { /* step through outedges of tail */
      if(node3!=head) CHANGE_STAT[0] += (OUT_DEG[node3] + IN_DEG[node3]);
     }
     STEP_THROUGH_INEDGES(tail, e, node3) { /* step through inedges of tail */
      if(node3!=head) CHANGE_STAT[0] += (OUT_DEG[node3] + IN_DEG[node3]);
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
// Rprintf("N_EDGES %d nedges %f \n",N_EDGES, nedges);
  CHANGE_STAT[0] /= nedges;
}

// A macro indicating whether x is in [from,to)
#define FROM_TO(x, from, to) ((x)>=(from) && (x)<(to))

/*****************
 changestat: d_degrange
*****************/
D_CHANGESTAT_FN(d_degrange) { 
  int i, j, echange;
  Vertex *id, *od;

  id=IN_DEG;
  od=OUT_DEG;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    Vertex tail, head;
    echange=(EdgetreeSearch(tail=TAIL(i), head=HEAD(i), nwp->outedges)==0)? 1:-1;
    Vertex taildeg = od[tail] + id[tail], headdeg = od[head] + id[head];
    for(j = 0; j < N_CHANGE_STATS; j++) {
      Vertex from = INPUT_PARAM[2*j], to = INPUT_PARAM[2*j+1];
      CHANGE_STAT[j] += FROM_TO(taildeg + echange, from, to) - FROM_TO(taildeg, from, to);
      CHANGE_STAT[j] += FROM_TO(headdeg + echange, from, to) - FROM_TO(headdeg, from, to);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}
 
/*****************
 changestat: d_degrange_by_attr
*****************/
D_CHANGESTAT_FN(d_degrange_by_attr) { 
  /* The inputparams are assumed to be set up as follows:
  The first 3*nstats values are in triples:  (from, to, attrvalue)
  The values following the first 3*nstats values are the nodal attributes.
  */
  int i, j;
  Vertex *id, *od;
  
  id=IN_DEG;
  od=OUT_DEG;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    Vertex tail, head;
    int echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1:1;
    Vertex taildeg = od[tail] + id[tail], headdeg = od[head] + id[head];
    int tailattr = INPUT_PARAM[3*N_CHANGE_STATS + tail - 1],
      headattr = INPUT_PARAM[3*N_CHANGE_STATS + head - 1]; 
    for(j = 0; j < N_CHANGE_STATS; j++) {
      Vertex from = INPUT_PARAM[3*j], to = INPUT_PARAM[3*j + 1];
      int testattr = INPUT_PARAM[3*j + 2]; 
      if (tailattr == testattr)  /* we have tail attr match */
        CHANGE_STAT[j] += FROM_TO(taildeg + echange, from, to) - FROM_TO(taildeg, from, to);
      if (headattr == testattr)  /* we have head attr match */
        CHANGE_STAT[j] += FROM_TO(headdeg + echange, from, to) - FROM_TO(headdeg, from, to);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_degrange_w_homophily
*****************/
D_CHANGESTAT_FN(d_degrange_w_homophily) { 
  /*  The inputparams are assumed to be set up as follows:
  The first 2*nstats values are the values of degrange
  The values following the first 2*nstats values are the nodal attributes.
  */
  int i, j;
  Vertex taildeg, headdeg, v;
  double *nodeattr;
  Edge e;

  nodeattr = mtp->inputparams + N_CHANGE_STATS*2 - 1;  

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    Vertex tail=TAIL(i), head=HEAD(i);
    int tailattr = nodeattr[tail], headattr = nodeattr[head];
    if (tailattr == headattr) { /* They match; otherwise don't bother */
      int echange = IS_OUTEDGE(tail, head) ? -1:1;
      taildeg=headdeg=-1; /* since tailattr==headattr, subtract the automatic match */
      taildeg=headdeg=0;
      STEP_THROUGH_OUTEDGES(tail, e, v) { taildeg += (nodeattr[v]==tailattr); }
      STEP_THROUGH_INEDGES(tail, e, v) { taildeg += (nodeattr[v]==tailattr); }
      STEP_THROUGH_OUTEDGES(head, e, v) { headdeg += (nodeattr[v]==headattr); }
      STEP_THROUGH_INEDGES(head, e, v) { headdeg += (nodeattr[v]==headattr); }
      for(j = 0; j < N_CHANGE_STATS; j++) {
        Vertex from = INPUT_PARAM[2*j], to = INPUT_PARAM[2*j+1];
        CHANGE_STAT[j] += FROM_TO(taildeg + echange, from, to) - FROM_TO(taildeg, from, to);
        CHANGE_STAT[j] += FROM_TO(headdeg + echange, from, to) - FROM_TO(headdeg, from, to);
      }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}                                        

#undef FROM_TO

/*****************
 changestat: d_degree
*****************/
D_CHANGESTAT_FN(d_degree) { 
  int i, j, echange;
  Vertex tail, head, taildeg, headdeg, deg, *id, *od;

  id=IN_DEG;
  od=OUT_DEG;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    echange=(EdgetreeSearch(tail=TAIL(i), head=HEAD(i), nwp->outedges)==0)? 1:-1;
    taildeg = od[tail] + id[tail];
    headdeg = od[head] + id[head];
    for(j = 0; j < N_CHANGE_STATS; j++) {
      deg = (Vertex)INPUT_PARAM[j];
      CHANGE_STAT[j] += (taildeg + echange == deg) - (taildeg == deg);
      CHANGE_STAT[j] += (headdeg + echange == deg) - (headdeg == deg);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}
 
/*****************
 changestat: d_degreepopularity
*****************/
D_CHANGESTAT_FN(d_degreepopularity) { 
  int i, edgeflag;
  double change;
  Vertex head, tail;
  
  /* *** don't forget tail -> head */    
  change = 0.0;
  FOR_EACH_TOGGLE(i) {
    tail=TAIL(i);
    head=HEAD(i);
    edgeflag = IS_UNDIRECTED_EDGE(tail, head); /* either 0 or 1 */
    Vertex tdeg = OUT_DEG[tail] + IN_DEG[tail];
    Vertex hdeg = OUT_DEG[head] + IN_DEG[head];
    if(edgeflag){
      change -= sqrt(tdeg);
      change += (tdeg-1.0)*(sqrt(tdeg-1.0)-sqrt(tdeg));
      change -= sqrt(hdeg);
      change += (hdeg-1.0)*(sqrt(hdeg-1.0)-sqrt(hdeg));
    }else{
      change += sqrt(tdeg+1.0);
      change += tdeg*(sqrt(tdeg+1.0)-sqrt(tdeg));
      change += sqrt(hdeg+1.0);
      change += hdeg*(sqrt(hdeg+1.0)-sqrt(hdeg));
    }
    TOGGLE_IF_MORE_TO_COME(i); 
  }
  CHANGE_STAT[0]=change; 
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_degree_by_attr
*****************/
D_CHANGESTAT_FN(d_degree_by_attr) { 
  /* The inputparams are assumed to be set up as follows:
  The first 2*nstats values are in pairs:  (degree, attrvalue)
  The values following the first 2*nstats values are the nodal attributes.
  */
  int i, j, echange, tailattr, headattr, testattr;
  Vertex tail, head, taildeg, headdeg, d, *id, *od;
  
  id=IN_DEG;
  od=OUT_DEG;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1:1;
    taildeg = od[tail] + id[tail];
    headdeg = od[head] + id[head];
    tailattr = INPUT_PARAM[2*N_CHANGE_STATS + tail - 1]; 
    headattr = INPUT_PARAM[2*N_CHANGE_STATS + head - 1]; 
    for(j = 0; j < N_CHANGE_STATS; j++) {
      d = (Vertex)INPUT_PARAM[2*j];
      testattr = INPUT_PARAM[2*j + 1]; 
      if (tailattr == testattr)  /* we have tail attr match */
        CHANGE_STAT[j] += (taildeg + echange == d) - (taildeg == d);
      if (headattr == testattr)  /* we have head attr match */
        CHANGE_STAT[j] += (headdeg + echange == d) - (headdeg == d);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_degree_w_homophily
*****************/
D_CHANGESTAT_FN(d_degree_w_homophily) { 
  /*  The inputparams are assumed to be set up as follows:
  The first nstats values are the values of degree
  The values following the first nstats values are the nodal attributes.
  */
  int i, j, echange, tailattr, headattr;
  Vertex tail, head, taildeg, headdeg, deg, v;
  double *nodeattr;
  Edge e;

  nodeattr = mtp->inputparams + N_CHANGE_STATS - 1;  

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    tail=TAIL(i);
    head=HEAD(i);
    tailattr = (int)nodeattr[tail];
    headattr = (int)nodeattr[head];    
    if (tailattr == headattr) { /* They match; otherwise don't bother */
      echange = IS_OUTEDGE(tail, head) ? -1:1;
      taildeg=headdeg=-1; /* since tailattr==headattr, subtract the automatic match */
      taildeg=headdeg=0;
      STEP_THROUGH_OUTEDGES(tail, e, v) { taildeg += (nodeattr[v]==tailattr); }
      STEP_THROUGH_INEDGES(tail, e, v) { taildeg += (nodeattr[v]==tailattr); }
      STEP_THROUGH_OUTEDGES(head, e, v) { headdeg += (nodeattr[v]==headattr); }
      STEP_THROUGH_INEDGES(head, e, v) { headdeg += (nodeattr[v]==headattr); }
      for(j = 0; j < N_CHANGE_STATS; j++) {
        deg = (Vertex)INPUT_PARAM[j];
        CHANGE_STAT[j] += (taildeg + echange == deg) - (taildeg == deg);
        CHANGE_STAT[j] += (headdeg + echange == deg) - (headdeg == deg);
      }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}                                        

/*****************
 changestat: d_density
*****************/
D_CHANGESTAT_FN(d_density) {
  int i;
  Edge ndyads = N_DYADS;
  
  /* *** don't forget tail -> head */    
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i) {
    CHANGE_STAT[0] += IS_OUTEDGE(TAIL(i), HEAD(i)) ? - 1 : 1;
    TOGGLE_IF_MORE_TO_COME(i);
  }
  CHANGE_STAT[0] = CHANGE_STAT[0] / ndyads;
  UNDO_PREVIOUS_TOGGLES(i);  
}

/*****************
 changestat: d_dsp
*****************/
D_CHANGESTAT_FN(d_dsp) { 
  Edge e, f;
  int i, j, echange;
  int L2tu, L2uh;
  Vertex deg;
  Vertex tail, head, u, v;
  
  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    echange = (EdgetreeSearch(tail=TAIL(i), head=HEAD(i), nwp->outedges) == 0) ? 1 : -1;
    /* step through outedges of head */
    for(e = EdgetreeMinimum(nwp->outedges, head);
    (u = nwp->outedges[e].value) != 0;
    e = EdgetreeSuccessor(nwp->outedges, e)){
      if (u != tail){
        L2tu=0;
        /* step through outedges of u */
        for(f = EdgetreeMinimum(nwp->outedges, u);
        (v = nwp->outedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->outedges, f)){
          if(EdgetreeSearch(MIN(v,tail),MAX(v,tail),nwp->outedges)!= 0) L2tu++;
        }
        /* step through inedges of u */
        for(f = EdgetreeMinimum(nwp->inedges, u); 
        (v = nwp->inedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->inedges, f)){
          if(EdgetreeSearch(MIN(v,tail),MAX(v,tail),nwp->outedges)!= 0) L2tu++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)INPUT_PARAM[j];
          CHANGE_STAT[j] += ((L2tu + echange == deg)
          - (L2tu == deg));
        }
      }
    }
    /* step through inedges of head */
    for(e = EdgetreeMinimum(nwp->inedges, head);
    (u = nwp->inedges[e].value) != 0;
    e = EdgetreeSuccessor(nwp->inedges, e)){
      if (u != tail){
        L2tu=0;
        /* step through outedges of u */
        for(f = EdgetreeMinimum(nwp->outedges, u);
        (v = nwp->outedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->outedges, f)){
          if(EdgetreeSearch(MIN(v,tail),MAX(v,tail),nwp->outedges)!= 0) L2tu++;
        }
        /* step through inedges of u */
        for(f = EdgetreeMinimum(nwp->inedges, u); 
        (v = nwp->inedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->inedges, f)){
          if(EdgetreeSearch(MIN(v,tail),MAX(v,tail),nwp->outedges)!= 0) L2tu++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)INPUT_PARAM[j];
          CHANGE_STAT[j] += ((L2tu + echange == deg)
          - (L2tu == deg));
        }
      }
    }
    /* step through outedges of tail */
    for(e = EdgetreeMinimum(nwp->outedges, tail);
    (u = nwp->outedges[e].value) != 0;
    e = EdgetreeSuccessor(nwp->outedges, e)){
      if (u != head){
        L2uh=0;
        /* step through outedges of u */
        for(f = EdgetreeMinimum(nwp->outedges, u);
        (v = nwp->outedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->outedges, f)){
          if(EdgetreeSearch(MIN(v,head),MAX(v,head),nwp->outedges)!= 0) L2uh++;
        }
        /* step through inedges of u */
        for(f = EdgetreeMinimum(nwp->inedges, u); 
        (v = nwp->inedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->inedges, f)){
          if(EdgetreeSearch(MIN(v,head),MAX(v,head),nwp->outedges)!= 0) L2uh++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)INPUT_PARAM[j];
          CHANGE_STAT[j] += ((L2uh + echange == deg)
          - (L2uh == deg));
        }
      }
    }
    /* step through inedges of tail */
    for(e = EdgetreeMinimum(nwp->inedges, tail);
    (u = nwp->inedges[e].value) != 0;
    e = EdgetreeSuccessor(nwp->inedges, e)){
      if (u != head){
        L2uh=0;
        /* step through outedges of u */
        for(f = EdgetreeMinimum(nwp->outedges, u);
        (v = nwp->outedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->outedges, f)){
          if(EdgetreeSearch(MIN(v,head),MAX(v,head),nwp->outedges)!= 0) L2uh++;
        }
        /* step through inedges of u */
        for(f = EdgetreeMinimum(nwp->inedges, u); 
        (v = nwp->inedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->inedges, f)){
          if(EdgetreeSearch(MIN(v,head),MAX(v,head),nwp->outedges)!= 0) L2uh++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)INPUT_PARAM[j];
          CHANGE_STAT[j] += ((L2uh + echange == deg)
          - (L2uh == deg));
        }
      }
    }    
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_dyadcov
*****************/
D_CHANGESTAT_FN(d_dyadcov) { 
  double val;
  Vertex tail, head;
  int i, edgeflag, refedgeflag;
  long int nrow, noffset, index;
  
  noffset = BIPARTITE;
  if(noffset > 0){
   nrow = (N_NODES)-(long int)(INPUT_PARAM[0]);
  }else{
   nrow = (long int)(INPUT_PARAM[0]);
  }
  
/*  Rprintf("nrow %d noffset %d\n",nrow, noffset);
  Rprintf("attrib: ");
  for(i=0;i<1000;i++)
   Rprintf("%1.0f",INPUT_ATTRIB[i]);

  Rprintf("\n;"); */

  if(DIRECTED){
    /* directed version */
    
    for(i=0;i<3;i++)
      CHANGE_STAT[i] = 0.0;
    
    /* *** don't forget tail -> head */    
    FOR_EACH_TOGGLE(i) {
      /*Get the initial state of the edge and its reflection*/
      edgeflag=IS_OUTEDGE(tail = TAIL(i), head = HEAD(i));
      refedgeflag = (EdgetreeSearch(head, tail, nwp->outedges) != 0);
      
      /*Get the dyadic covariate*/
      /*    val = INPUT_ATTRIB[(head-1-nrow)+(tail-1)*ncols]; */
      index = (head-1-noffset)*nrow+(tail-1);
      if(index >= 0 && index <= nrow*nrow){
        val = INPUT_ATTRIB[(head-1-noffset)*nrow+(tail-1)];
        /*  Rprintf("tail %d head %d nrow %d ncols %d val %f\n",tail, head, nrow, ncols, val); */
        
        /*Update the change statistics, as appropriate*/
        if(refedgeflag){      /* Reflected edge is present */
          if(edgeflag){         /* Toggled edge _was_ present */
            if(head>tail){              /* Mut to low->high */
              CHANGE_STAT[0] -= val;
              CHANGE_STAT[1] += val;
            }else{                /* Mut to high->low */
              CHANGE_STAT[0] -= val;
              CHANGE_STAT[2] += val;
            }
          }else{                /* Toggled edge _was not_ present */
            if(head>tail){              /* Low->high to mut */
              CHANGE_STAT[1] -= val;
              CHANGE_STAT[0] += val;
            }else{                /* High->low to mut */
              CHANGE_STAT[2] -= val;
              CHANGE_STAT[0] += val;
            }
          }
        }else{                /* Reflected edge is absent */
          if(edgeflag){         /* Toggled edge _was_ present */
            if(head>tail){              /* High->low to null */
              CHANGE_STAT[2] -= val;
            }else{                /* Low->high to null */
              CHANGE_STAT[1] -= val;
            }
          }else{                /* Toggled edge _was not_ present */
            if(head>tail){              /* Null to high->low */
              CHANGE_STAT[2] += val;
            }else{                /* Null to low->high */
              CHANGE_STAT[1] += val;
            }
          }
        }
      }      
      TOGGLE_IF_MORE_TO_COME(i);
    }
  }else{
    /* undirected case (including bipartite) */
    
    /* *** don't forget tail -> head */    
    CHANGE_STAT[0] = 0.0;
    FOR_EACH_TOGGLE(i) {
      /*Get the initial edge state*/
      edgeflag=IS_OUTEDGE(tail = TAIL(i), head = HEAD(i));
      /*Get the covariate value*/
      /*    val = INPUT_ATTRIB[(head-1-nrow)+(tail-1)*ncols]; */
      index = (head-1-noffset)*nrow+(tail-1);
      if(index >= 0 && index <= nrow*((long int)(INPUT_PARAM[0]))){
        val = INPUT_ATTRIB[(head-1-noffset)*nrow+(tail-1)];
        /*Update the change statistic, based on the toggle type*/
        /*  Rprintf("tail %d head %d nrow %d noffset %d val %f\n",tail, head, nrow, noffset, val); */
        /*Update the change statistic, based on the toggle type*/
        CHANGE_STAT[0] += edgeflag ? -val : val;
      }
      TOGGLE_IF_MORE_TO_COME(i);
    }
  }
  UNDO_PREVIOUS_TOGGLES(i);
}


/********************  changestats:  E    ***********/
/*****************
 changestat: d_edgecov
*****************/
D_CHANGESTAT_FN(d_edgecov) {
  double val;
  Vertex tail, head;
  int nrow, noffset;
  int i, edgeflag;
  
  noffset = BIPARTITE;
  if(noffset > 0){
    /*   nrow = (N_NODES)-(long int)(INPUT_PARAM[0]); */
    nrow = noffset;
  }else{
    nrow = (long int)(INPUT_PARAM[0]);
  }
  
  /* *** don't forget tail -> head */    
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i) {
    /*Get the initial edge state*/
    edgeflag=IS_OUTEDGE(tail=TAIL(i), head=HEAD(i));
    /*Get the covariate value*/
    /*    val = INPUT_ATTRIB[(head-1-nrow)+(tail-1)*ncols]; */
    /*OLD COMMENT:  Note: h/t are backwards!*/
    /* *** despite the comment above, I swapped heads/tails in the line below. is this right? */
    val = INPUT_ATTRIB[(head-1-noffset)*nrow+(tail-1)];  
    /*  Rprintf("tail %d head %d nrow %d val %f\n", tail, head, nrow, val); */
    /*Update the change statistic, based on the toggle type*/
    CHANGE_STAT[0] += edgeflag ? -val : val;
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_edges
*****************/
D_CHANGESTAT_FN(d_edges) {
  int edgeflag, i;
  Vertex tail, head;
  
  /* *** don't forget tail -> head */    
  CHANGE_STAT[0] = 0.0;
  for (i=0; i < ntoggles; i++)
    {
      edgeflag = IS_OUTEDGE(tail = TAIL(i), head = HEAD(i));
      CHANGE_STAT[0] += edgeflag ? - 1 : 1;
      TOGGLE_IF_MORE_TO_COME(i);
    }
  UNDO_PREVIOUS_TOGGLES(i);
}

S_CHANGESTAT_FN(s_edges) {
  CHANGE_STAT[0] = N_EDGES;
}

/*****************
 changestat: d_esp
*****************/
D_CHANGESTAT_FN(d_esp) { 
  Edge e, f;
  int i, j, echange;
  int L2th, L2tu, L2uh;
  Vertex deg;
  Vertex tail, head, u, v;
  
  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    L2th=0;
    echange = (EdgetreeSearch(tail=TAIL(i), head=HEAD(i), nwp->outedges) == 0) ? 1 : -1;
    /* step through outedges of head */
    for(e = EdgetreeMinimum(nwp->outedges, head); 
    (u = nwp->outedges[e].value) != 0; e = EdgetreeSuccessor(nwp->outedges, e)) {
      if (EdgetreeSearch(MIN(u,tail), MAX(u,tail), nwp->outedges) != 0){
        L2th++;
        L2tu=0;
        L2uh=0;
        /* step through outedges of u */
        for(f = EdgetreeMinimum(nwp->outedges, u); (v = nwp->outedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->outedges, f)){
          if(EdgetreeSearch(MIN(v,head),MAX(v,head),nwp->outedges)!= 0) L2uh++;
          if(EdgetreeSearch(MIN(v,tail),MAX(v,tail),nwp->outedges)!= 0) L2tu++;
        }
        /* step through inedges of u */
        for(f = EdgetreeMinimum(nwp->inedges, u); (v = nwp->inedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->inedges, f)){
          if(EdgetreeSearch(MIN(v,head),MAX(v,head),nwp->outedges)!= 0) L2uh++;
          if(EdgetreeSearch(MIN(v,tail),MAX(v,tail),nwp->outedges)!= 0) L2tu++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)INPUT_PARAM[j];
          CHANGE_STAT[j] += ((L2tu + echange == deg)
          - (L2tu == deg));
          CHANGE_STAT[j] += ((L2uh + echange == deg)
          - (L2uh == deg));
        }
      }
    }
    /* step through inedges of head */
    for (e = EdgetreeMinimum(nwp->inedges, head); (u = nwp->inedges[e].value) != 0;
    e = EdgetreeSuccessor(nwp->inedges, e)){
      if (EdgetreeSearch(MIN(u,tail), MAX(u,tail), nwp->outedges) != 0){
        L2th++;
        L2tu=0;
        L2uh=0;
        /* step through outedges of u */
        for(f = EdgetreeMinimum(nwp->outedges, u); (v = nwp->outedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->outedges, f)){
          if(EdgetreeSearch(MIN(v,head),MAX(v,head),nwp->outedges)!= 0) L2uh++;
          if(EdgetreeSearch(MIN(v,tail),MAX(v,tail),nwp->outedges)!= 0) L2tu++;
        }
        /* step through inedges of u */
        for(f = EdgetreeMinimum(nwp->inedges, u); (v = nwp->inedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->inedges, f)){
          if(EdgetreeSearch(MIN(v,head),MAX(v,head),nwp->outedges)!= 0) L2uh++;
          if(EdgetreeSearch(MIN(v,tail),MAX(v,tail),nwp->outedges)!= 0) L2tu++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)INPUT_PARAM[j];
          CHANGE_STAT[j] += ((L2tu + echange == deg)
          - (L2tu == deg));
          CHANGE_STAT[j] += ((L2uh + echange == deg)
          - (L2uh == deg));
        }
      }
    }
    for(j = 0; j < N_CHANGE_STATS; j++){
      deg = (Vertex)INPUT_PARAM[j];
/*      CHANGE_STAT[j] += echange*((L2th == deg) - (0 == deg)); */
      CHANGE_STAT[j] += echange*(L2th == deg);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }  
  UNDO_PREVIOUS_TOGGLES(i);
}

/********************  changestats:  F    ***********/

/********************  changestats:  G    ***********/
/*****************
 changestat: d_gwb1degree
*****************/
D_CHANGESTAT_FN(d_gwb1degree) { 
  /* It is assumed that in this bipartite network, the only edges are
  of the form (b1, b2), where b1 is always strictly less
  than b2.  In other words, the degree of a b1 is equivalent
  to its outdegree and the degree of a b2 is equivalent to its
  indegree.
  */
  int i, echange;
  double decay, oneexpd;
  Vertex b1, b1deg, *od;
  
  decay = INPUT_PARAM[0];
  oneexpd = 1.0-exp(-decay);
  od=OUT_DEG;

  /* *** don't forget tail -> head */    
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i) {
    echange=(EdgetreeSearch(b1=TAIL(i), HEAD(i), nwp->outedges)==0) ? 1 : -1;
    b1deg = od[b1]+(echange-1)/2;
    CHANGE_STAT[0] += echange*pow(oneexpd,(double)b1deg);
    TOGGLE_IF_MORE_TO_COME(i);
  }
  
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_gwb1degree_by_attr
*****************/
D_CHANGESTAT_FN(d_gwb1degree_by_attr) { 
  /* It is assumed that in this bipartite network, the only edges are
  of the form (b1, b2), where b1 is always strictly less
  than b2.  In other words, the degree of a b1 is equivalent
  to its outdegree and the degree of a b2 is equivalent to its
  indegree.
  The inputparams are assumed to be set up as follows:
    The first value is theta (as in Hunter et al, JASA 200?), controlling decay
    The next sequence of values is the nodal attributes, coded as integers
         from 1 through N_CHANGE_STATS
  */
  int i, echange, b1attr;
  double decay, oneexpd;
  Vertex b1, b1deg, *od;
  
  decay = INPUT_PARAM[0];
  oneexpd = 1.0-exp(-decay);
  od=OUT_DEG;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    echange=(EdgetreeSearch(b1=TAIL(i), HEAD(i), nwp->outedges)==0) ? 1 : -1;
    b1deg = od[b1]+(echange-1)/2;
    b1attr = INPUT_PARAM[b1]; 
    /* *** the comment below looked right, so I didn't swap it - ALC */
    /*  Rprintf("b1 %d heads %d b1deg %d b1attr %d echange %d\n",b1, HEAD(i), b1deg, b1attr, echange); */
    CHANGE_STAT[b1attr-1] += echange * pow(oneexpd,(double)b1deg);
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_gwdegree
*****************/
D_CHANGESTAT_FN(d_gwdegree) { 
  int i, echange=0;
  double decay, oneexpd, change;
  Vertex tail, head, taild, headd=0, *id, *od;
  
  id=IN_DEG;
  od=OUT_DEG;
  decay = INPUT_PARAM[0];
  oneexpd = 1.0-exp(-decay);
  
  /* *** don't forget tail -> head */    
  change = 0.0;
  FOR_EACH_TOGGLE(i) {      
    echange = (EdgetreeSearch(tail=TAIL(i), head=HEAD(i), nwp->outedges) == 0) ? 1 : -1;
    taild = od[tail] + id[tail] + (echange - 1)/2;
    headd = od[head] + id[head] + (echange - 1)/2;
    change += echange*(pow(oneexpd,(double)taild)+pow(oneexpd,(double)headd));
      
    TOGGLE_IF_MORE_TO_COME(i);
  }
  CHANGE_STAT[0] = change;
  
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_gwdegree_by_attr
*****************/
D_CHANGESTAT_FN(d_gwdegree_by_attr) { 
  /*The inputparams are assumed to be set up as follows:
    The first value is the decay parameter (as in Hunter et al, JASA 200?)
    The next sequence of values is the nodal attributes, coded as integers
         from 1 through N_CHANGE_STATS
  */
  int i, tailattr, headattr, echange=0;
  double decay, oneexpd;
  Vertex tail, head, taild, headd=0, *id, *od;
  
  id=IN_DEG;
  od=OUT_DEG;
  decay = INPUT_PARAM[0];
  oneexpd = 1.0-exp(-decay);
  
  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    echange = (EdgetreeSearch(tail=TAIL(i), head=HEAD(i), nwp->outedges) == 0) ? 1 : -1;
    taild = od[tail] + id[tail] + (echange - 1)/2;
    tailattr = INPUT_PARAM[tail]; 
    CHANGE_STAT[tailattr-1] += echange*(pow(oneexpd,(double)taild));
    
    headd = od[head] + id[head] + (echange - 1)/2;
    headattr = INPUT_PARAM[head]; 
    CHANGE_STAT[headattr-1] += echange*(pow(oneexpd,(double)headd));
      
    TOGGLE_IF_MORE_TO_COME(i);
  }
  
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_gwdsp
****************/
D_CHANGESTAT_FN(d_gwdsp) { 
  Edge e, f;
  int i, echange, ochange;
  int L2tu, L2uh;
  Vertex tail, head, u, v;
  double alpha, oneexpa, cumchange;
  
  /* *** don't forget tail -> head */    
  CHANGE_STAT[0] = 0.0;
  alpha = INPUT_PARAM[0];
  oneexpa = 1.0-exp(-alpha);
  
  FOR_EACH_TOGGLE(i){      
    cumchange=0.0;
    ochange = (EdgetreeSearch(tail=TAIL(i), head=HEAD(i), nwp->outedges) == 0) ? 0 : -1;
    echange = 2*ochange + 1;
    /* step through outedges of head */
    for(e = EdgetreeMinimum(nwp->outedges, head);
    (u = nwp->outedges[e].value) != 0;
    e = EdgetreeSuccessor(nwp->outedges, e)){
      if (u != tail){
        L2tu=ochange;
        /* step through outedges of u */
        for(f = EdgetreeMinimum(nwp->outedges, u);
        (v = nwp->outedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->outedges, f)){
          if(EdgetreeSearch(MIN(v,tail),MAX(v,tail),nwp->outedges)!= 0) L2tu++;
        }
        /* step through inedges of u */
        for(f = EdgetreeMinimum(nwp->inedges, u);
        (v = nwp->inedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->inedges, f)){
          if(EdgetreeSearch(MIN(v,tail),MAX(v,tail),nwp->outedges)!= 0) L2tu++;
        }
        cumchange += pow(oneexpa,(double)L2tu);
      }
    }
    /* step through inedges of head */
    for(e = EdgetreeMinimum(nwp->inedges, head);
    (u = nwp->inedges[e].value) != 0;
    e = EdgetreeSuccessor(nwp->inedges, e)){
      if (u != tail){
        L2tu=ochange;
        /* step through outedges of u */
        for(f = EdgetreeMinimum(nwp->outedges, u);
        (v = nwp->outedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->outedges, f)){
          if(EdgetreeSearch(MIN(v,tail),MAX(v,tail),nwp->outedges)!= 0) L2tu++;
        }
        /* step through inedges of u */
        for(f = EdgetreeMinimum(nwp->inedges, u);
        (v = nwp->inedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->inedges, f)){
          if(EdgetreeSearch(MIN(v,tail),MAX(v,tail),nwp->outedges)!= 0) L2tu++;
        }
        cumchange += pow(oneexpa,(double)L2tu);
      }
    }
    
    /* step through outedges of tail  */
    for(e = EdgetreeMinimum(nwp->outedges, tail);
    (u = nwp->outedges[e].value) != 0;
    e = EdgetreeSuccessor(nwp->outedges, e)){
      if (u != head){
        L2uh=ochange;
        /* step through outedges of u */
        for(f = EdgetreeMinimum(nwp->outedges, u);
        (v = nwp->outedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->outedges, f)){
          if(EdgetreeSearch(MIN(v,head),MAX(v,head),nwp->outedges)!= 0) L2uh++;
        }
        /* step through inedges of u */
        for(f = EdgetreeMinimum(nwp->inedges, u); 
        (v = nwp->inedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->inedges, f)){
          if(EdgetreeSearch(MIN(v,head),MAX(v,head),nwp->outedges)!= 0) L2uh++;
        }
        cumchange += pow(oneexpa,(double)L2uh);
      }
    }
    /* step through inedges of tail */
    for(e = EdgetreeMinimum(nwp->inedges, tail);
    (u = nwp->inedges[e].value) != 0;
    e = EdgetreeSuccessor(nwp->inedges, e)){
      if (u != head){
        L2uh=ochange;
        /* step through outedges of u */
        for(f = EdgetreeMinimum(nwp->outedges, u);
        (v = nwp->outedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->outedges, f)){
          if(EdgetreeSearch(MIN(v,head),MAX(v,head),nwp->outedges)!= 0) L2uh++;
        }
        /* step through inedges of u */
        for(f = EdgetreeMinimum(nwp->inedges, u); 
        (v = nwp->inedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->inedges, f)){
          if(EdgetreeSearch(MIN(v,head),MAX(v,head),nwp->outedges)!= 0) L2uh++;
        }
        cumchange += pow(oneexpa,(double)L2uh);
      }
    }
    
    cumchange  = echange*cumchange;
    (CHANGE_STAT[0]) += cumchange;
    TOGGLE_IF_MORE_TO_COME(i);
  }
  
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_gwb2degree
*****************/
D_CHANGESTAT_FN(d_gwb2degree) { 
  /* It is assumed that in this bipartite network, the only edges are
  of the form (b1, b2), where b1 is always strictly less
  than b2.  In other words, the degree of a b1 is equivalent
  to its outdegree and the degree of a b2 is equivalent to its
  indegree.
  */
  int i, echange;
  double decay, oneexpd;
  Vertex b2, b2deg, *id;
  
  decay = INPUT_PARAM[0];
  oneexpd = 1.0-exp(-decay);
  id=IN_DEG;

  /* *** don't forget tail -> head */    
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i) {      
    echange=(EdgetreeSearch(TAIL(i), b2=HEAD(i), nwp->outedges)==0) ? 1 : -1;
    b2deg = id[b2]+(echange-1)/2;
    CHANGE_STAT[0] += echange*pow(oneexpd,(double)b2deg);
    TOGGLE_IF_MORE_TO_COME(i);
  }
  
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_gwb2degree_by_attr
*****************/
D_CHANGESTAT_FN(d_gwb2degree_by_attr) { 
  /* It is assumed that in this bipartite network, the only edges are
  of the form (b1, b2), where b1 is always strictly less
  than b2.  In other words, the degree of a b1 is equivalent
  to its outdegree and the degree of a b2 is equivalent to its
  indegree.
  The inputparams are assumed to be set up as follows:
    The first value is theta (as in Hunter et al, JASA 200?), controlling decay
    The next sequence of values is the nodal attributes, coded as integers
         from 1 through N_CHANGE_STATS
  */
  int i, echange, b2attr;
  double decay, oneexpd;
  Vertex b2, b2deg, *id;
  
  decay = INPUT_PARAM[0];
  oneexpd = 1.0-exp(-decay);
  id=IN_DEG;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {      
    echange=(EdgetreeSearch(TAIL(i), b2=HEAD(i), nwp->outedges)==0) ? 1 : -1;
    b2deg = id[b2]+(echange-1)/2;
    b2attr = INPUT_PARAM[b2]; 
/*  Rprintf("tail %d b2 %d b2deg %d b2attr %d echange %d\n",TAIL(i), b2, b2deg, b2attr, echange); */
    CHANGE_STAT[b2attr-1] += echange * pow(oneexpd,(double)b2deg);
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_gwesp
*****************/
D_CHANGESTAT_FN(d_gwesp) { 
  Edge e, f;
  int i, echange, ochange;
  int L2th, L2tu, L2uh;
  Vertex tail, head, u, v;
  double alpha, oneexpa, cumchange;
  
  CHANGE_STAT[0] = 0.0;
  alpha = INPUT_PARAM[0];
  oneexpa = 1.0-exp(-alpha);
  
  /* *** don't forget tail -> head */    
  FOR_EACH_TOGGLE(i){
    cumchange=0.0;
    L2th=0;
    ochange = (EdgetreeSearch(tail=TAIL(i), head=HEAD(i), nwp->outedges) == 0) ? 0 : -1;
    echange = 2*ochange + 1;
    /* step through outedges of head  */
    for(e = EdgetreeMinimum(nwp->outedges, head);
	(u = nwp->outedges[e].value) != 0;
	e = EdgetreeSuccessor(nwp->outedges, e)){
      if (EdgetreeSearch(MIN(u,tail), MAX(u,tail), nwp->outedges) != 0){
	L2th++;
	L2tu=ochange;
	L2uh=ochange;
	/* step through outedges of u */
	for(f = EdgetreeMinimum(nwp->outedges, u);
	    (v = nwp->outedges[f].value) != 0;
	    f = EdgetreeSuccessor(nwp->outedges, f)){
	  if(EdgetreeSearch(MIN(v,head),MAX(v,head),nwp->outedges)!= 0) L2uh++;
	  if(EdgetreeSearch(MIN(v,tail),MAX(v,tail),nwp->outedges)!= 0) L2tu++;
	}
	/* step through inedges of u */
	for(f = EdgetreeMinimum(nwp->inedges, u); 
	    (v = nwp->inedges[f].value) != 0;
	    f = EdgetreeSuccessor(nwp->inedges, f)){
	  if(EdgetreeSearch(MIN(v,head),MAX(v,head),nwp->outedges)!= 0) L2uh++;
	  if(EdgetreeSearch(MIN(v,tail),MAX(v,tail),nwp->outedges)!= 0) L2tu++;
	}
	cumchange += pow(oneexpa,(double)L2tu) +
	  pow(oneexpa,(double)L2uh) ;
      }
    }
    /* step through inedges of head */
    
    for(e = EdgetreeMinimum(nwp->inedges, head);
	(u = nwp->inedges[e].value) != 0;
	e = EdgetreeSuccessor(nwp->inedges, e)){
      if (EdgetreeSearch(MIN(u,tail), MAX(u,tail), nwp->outedges) != 0){
	L2th++;
	L2tu=ochange;
	L2uh=ochange;
	/* step through outedges of u */
	for(f = EdgetreeMinimum(nwp->outedges, u);
	    (v = nwp->outedges[f].value) != 0;
	    f = EdgetreeSuccessor(nwp->outedges, f)){
	  if(EdgetreeSearch(MIN(v,head),MAX(v,head),nwp->outedges)!= 0) L2uh++;
	  if(EdgetreeSearch(MIN(v,tail),MAX(v,tail),nwp->outedges)!= 0) L2tu++;
	}
	/* step through inedges of u */
	for(f = EdgetreeMinimum(nwp->inedges, u); 
	    (v = nwp->inedges[f].value) != 0;
	    f = EdgetreeSuccessor(nwp->inedges, f)){
	  if(EdgetreeSearch(MIN(v,head),MAX(v,head),nwp->outedges)!= 0) L2uh++;
	  if(EdgetreeSearch(MIN(v,tail),MAX(v,tail),nwp->outedges)!= 0) L2tu++;
	}
	cumchange += pow(oneexpa,(double)L2tu) +
	  pow(oneexpa,(double)L2uh) ;
      }
    }
    
    if(alpha < 100.0){
      cumchange += exp(alpha)*(1.0-pow(oneexpa,(double)L2th)) ;
    }else{
      cumchange += (double)L2th;
    }
    cumchange  = echange*cumchange;
    (CHANGE_STAT[0]) += cumchange;
    TOGGLE_IF_MORE_TO_COME(i);
  }
  
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_gwidegree
*****************/
D_CHANGESTAT_FN(d_gwidegree) { 
  int i, edgeflag;
  double decay, oneexpd, change;
  Vertex head, headd=0;
  
  decay = INPUT_PARAM[0];
  oneexpd = 1.0-exp(-decay);
  change = 0.0;

  /* *** don't forget tail -> head */    
  FOR_EACH_TOGGLE(i) {
    head=HEAD(i);
    edgeflag = IS_OUTEDGE(TAIL(i), head); /* either 0 or 1 */
    headd = IN_DEG[head] - edgeflag;
    change += (edgeflag? -1.0 : 1.0) * pow(oneexpd,(double)headd);
    TOGGLE_IF_MORE_TO_COME(i); 
  }
  CHANGE_STAT[0]=change; 
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_gwidegree_by_attr
*****************/
D_CHANGESTAT_FN(d_gwidegree_by_attr) { 
  /*The inputparams are assumed to be set up as follows:
    The first value is the decay parameter (as in Hunter et al, JASA 2008)
    The next sequence of values is the nodal attributes, coded as integers
         from 1 through N_CHANGE_STATS
  */
  int i, headattr, echange;
  double decay, oneexpd;
  Vertex head, headd;
  
  decay = INPUT_PARAM[0];
  oneexpd = 1.0-exp(-decay);

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    head = HEAD(i);
    echange = IS_OUTEDGE(TAIL(i), head) ? -1 : 1;
    headd = IN_DEG[head] + (echange - 1)/2;
    headattr = INPUT_PARAM[head]; 
    CHANGE_STAT[headattr-1] += echange*(pow(oneexpd,(double)headd));      
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_gwnsp
*****************/
D_CHANGESTAT_FN(d_gwnsp) { 
  Edge e, f;
  int i, echange, ochange;
  int L2th, L2tu, L2uh;
  Vertex tail, head, u, v;
  double alpha, oneexpa, cumchange;
  
  CHANGE_STAT[0] = 0.0;

  alpha = INPUT_PARAM[0];
  oneexpa = 1.0-exp(-alpha);
  
  /* *** don't forget tail -> head */    
  FOR_EACH_TOGGLE(i){      
    cumchange=0.0;
    ochange = (EdgetreeSearch(tail=TAIL(i), head=HEAD(i), nwp->outedges) == 0) ? 0 : -1;
    echange = 2*ochange + 1;
    /* step through outedges of head */
    for(e = EdgetreeMinimum(nwp->outedges, head);
    (u = nwp->outedges[e].value) != 0;
    e = EdgetreeSuccessor(nwp->outedges, e)){
      if (u != tail){
        L2tu=ochange;
        /* step through outedges of u */
        for(f = EdgetreeMinimum(nwp->outedges, u);
        (v = nwp->outedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->outedges, f)){
          if(EdgetreeSearch(MIN(v,tail),MAX(v,tail),nwp->outedges)!= 0) L2tu++;
        }
        /* step through inedges of u */
        for(f = EdgetreeMinimum(nwp->inedges, u);
        (v = nwp->inedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->inedges, f)){
          if(EdgetreeSearch(MIN(v,tail),MAX(v,tail),nwp->outedges)!= 0) L2tu++;
        }
        cumchange += pow(oneexpa,(double)L2tu);
      }
    }
    /* step through inedges of head */
    for(e = EdgetreeMinimum(nwp->inedges, head);
    (u = nwp->inedges[e].value) != 0;
    e = EdgetreeSuccessor(nwp->inedges, e)){
      if (u != tail){
        L2tu=ochange;
        /* step through outedges of u */
        for(f = EdgetreeMinimum(nwp->outedges, u);
        (v = nwp->outedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->outedges, f)){
          if(EdgetreeSearch(MIN(v,tail),MAX(v,tail),nwp->outedges)!= 0) L2tu++;
        }
        /* step through inedges of u */
        for(f = EdgetreeMinimum(nwp->inedges, u);
        (v = nwp->inedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->inedges, f)){
          if(EdgetreeSearch(MIN(v,tail),MAX(v,tail),nwp->outedges)!= 0) L2tu++;
        }
        cumchange += pow(oneexpa,(double)L2tu);
      }
    }
    
    /* step through outedges of tail  */
    for(e = EdgetreeMinimum(nwp->outedges, tail);
    (u = nwp->outedges[e].value) != 0;
    e = EdgetreeSuccessor(nwp->outedges, e)){
      if (u != head){
        L2uh=ochange;
        /* step through outedges of u */
        for(f = EdgetreeMinimum(nwp->outedges, u);
        (v = nwp->outedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->outedges, f)){
          if(EdgetreeSearch(MIN(v,head),MAX(v,head),nwp->outedges)!= 0) L2uh++;
        }
        /* step through inedges of u */
        for(f = EdgetreeMinimum(nwp->inedges, u); 
        (v = nwp->inedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->inedges, f)){
          if(EdgetreeSearch(MIN(v,head),MAX(v,head),nwp->outedges)!= 0) L2uh++;
        }
        cumchange += pow(oneexpa,(double)L2uh);
      }
    }
    /* step through inedges of tail */
    for(e = EdgetreeMinimum(nwp->inedges, tail);
    (u = nwp->inedges[e].value) != 0;
    e = EdgetreeSuccessor(nwp->inedges, e)){
      if (u != head){
        L2uh=ochange;
        /* step through outedges of u */
        for(f = EdgetreeMinimum(nwp->outedges, u);
        (v = nwp->outedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->outedges, f)){
          if(EdgetreeSearch(MIN(v,head),MAX(v,head),nwp->outedges)!= 0) L2uh++;
        }
        /* step through inedges of u */
        for(f = EdgetreeMinimum(nwp->inedges, u); 
        (v = nwp->inedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->inedges, f)){
          if(EdgetreeSearch(MIN(v,head),MAX(v,head),nwp->outedges)!= 0) L2uh++;
        }
        cumchange += pow(oneexpa,(double)L2uh);
      }
    }
    
    cumchange  = echange*cumchange;
    (CHANGE_STAT[0]) += cumchange;
    TOGGLE_IF_MORE_TO_COME(i);
  }
  
  UNDO_PREVIOUS_TOGGLES(i);

  alpha = INPUT_PARAM[0];
  oneexpa = 1.0-exp(-alpha);
  
  /* *** don't forget tail -> head */    
  FOR_EACH_TOGGLE(i){
    cumchange=0.0;
    L2th=0;
    ochange = (EdgetreeSearch(tail=TAIL(i), head=HEAD(i), nwp->outedges) == 0) ? 0 : -1;
    echange = 2*ochange + 1;
    /* step through outedges of head  */
    for(e = EdgetreeMinimum(nwp->outedges, head);
	(u = nwp->outedges[e].value) != 0;
	e = EdgetreeSuccessor(nwp->outedges, e)){
      if (EdgetreeSearch(MIN(u,tail), MAX(u,tail), nwp->outedges) != 0){
	L2th++;
	L2tu=ochange;
	L2uh=ochange;
	/* step through outedges of u */
	for(f = EdgetreeMinimum(nwp->outedges, u);
	    (v = nwp->outedges[f].value) != 0;
	    f = EdgetreeSuccessor(nwp->outedges, f)){
	  if(EdgetreeSearch(MIN(v,head),MAX(v,head),nwp->outedges)!= 0) L2uh++;
	  if(EdgetreeSearch(MIN(v,tail),MAX(v,tail),nwp->outedges)!= 0) L2tu++;
	}
	/* step through inedges of u */
	for(f = EdgetreeMinimum(nwp->inedges, u); 
	    (v = nwp->inedges[f].value) != 0;
	    f = EdgetreeSuccessor(nwp->inedges, f)){
	  if(EdgetreeSearch(MIN(v,head),MAX(v,head),nwp->outedges)!= 0) L2uh++;
	  if(EdgetreeSearch(MIN(v,tail),MAX(v,tail),nwp->outedges)!= 0) L2tu++;
	}
	cumchange += pow(oneexpa,(double)L2tu) +
	  pow(oneexpa,(double)L2uh) ;
      }
    }
    /* step through inedges of head */
    
    for(e = EdgetreeMinimum(nwp->inedges, head);
	(u = nwp->inedges[e].value) != 0;
	e = EdgetreeSuccessor(nwp->inedges, e)){
      if (EdgetreeSearch(MIN(u,tail), MAX(u,tail), nwp->outedges) != 0){
	L2th++;
	L2tu=ochange;
	L2uh=ochange;
	/* step through outedges of u */
	for(f = EdgetreeMinimum(nwp->outedges, u);
	    (v = nwp->outedges[f].value) != 0;
	    f = EdgetreeSuccessor(nwp->outedges, f)){
	  if(EdgetreeSearch(MIN(v,head),MAX(v,head),nwp->outedges)!= 0) L2uh++;
	  if(EdgetreeSearch(MIN(v,tail),MAX(v,tail),nwp->outedges)!= 0) L2tu++;
	}
	/* step through inedges of u */
	for(f = EdgetreeMinimum(nwp->inedges, u); 
	    (v = nwp->inedges[f].value) != 0;
	    f = EdgetreeSuccessor(nwp->inedges, f)){
	  if(EdgetreeSearch(MIN(v,head),MAX(v,head),nwp->outedges)!= 0) L2uh++;
	  if(EdgetreeSearch(MIN(v,tail),MAX(v,tail),nwp->outedges)!= 0) L2tu++;
	}
	cumchange += pow(oneexpa,(double)L2tu) +
	  pow(oneexpa,(double)L2uh) ;
      }
    }
    
    if(alpha < 100.0){
      cumchange += exp(alpha)*(1.0-pow(oneexpa,(double)L2th)) ;
    }else{
      cumchange += (double)L2th;
    }
    cumchange  = echange*cumchange;
    (CHANGE_STAT[0]) -= cumchange;
    TOGGLE_IF_MORE_TO_COME(i);
  }
  
  UNDO_PREVIOUS_TOGGLES(i);

}

/*****************
 changestat: d_gwodegree
*****************/
D_CHANGESTAT_FN(d_gwodegree) { 
  int i, edgeflag;
  double decay, oneexpd, change;
  Vertex tail, taild;
  
  decay = INPUT_PARAM[0];
  oneexpd = 1.0-exp(-decay);  
  change = 0.0;

  /* *** don't forget tail -> head */    
  FOR_EACH_TOGGLE(i) {
    tail=TAIL(i);
    edgeflag = IS_OUTEDGE(tail, HEAD(i)); /* either 0 or 1 */
    taild = OUT_DEG[tail] - edgeflag;
    change += (edgeflag? -1 : 1) * pow(oneexpd,(double)taild);
    TOGGLE_IF_MORE_TO_COME(i);
  }
  CHANGE_STAT[0] = change;
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_gwodegree_by_attr
*****************/
D_CHANGESTAT_FN(d_gwodegree_by_attr) { 
  /*The inputparams are assumed to be set up as follows:
    The first value is the decay parameter (as in Hunter et al, JASA 2008)
    The next sequence of values is the nodal attributes, coded as integers
         from 1 through N_CHANGE_STATS
  */
  int i, tailattr, echange;
  double decay, oneexpd;
  Vertex tail, taild;
  
  decay = INPUT_PARAM[0];
  oneexpd = 1.0-exp(-decay);

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    tail = TAIL(i);
    echange = IS_OUTEDGE(tail, HEAD(i)) ? -1 : 1;
    taild = OUT_DEG[tail] + (echange - 1)/2;
    tailattr = INPUT_PARAM[tail]; 
    CHANGE_STAT[tailattr-1] += echange*(pow(oneexpd,(double)taild));      
    TOGGLE_IF_MORE_TO_COME(i);
  }  
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_gwtdsp
****************/
D_CHANGESTAT_FN(d_gwtdsp) {
  Edge e, f;
  int i, echange, ochange, L2tu, L2uh;
  Vertex tail, head, u, v;
  double alpha, oneexpa, cumchange;
  
  CHANGE_STAT[0] = 0.0;
  alpha = INPUT_PARAM[0];
  oneexpa = 1.0-exp(-alpha);
  
  /* *** don't forget tail -> head */    
  FOR_EACH_TOGGLE(i){
    tail=TAIL(i); head=HEAD(i);
    cumchange=0.0;
    ochange = -IS_OUTEDGE(tail,head);
    echange = 2*ochange + 1;
    /* step through outedges of head */
    for(e = MIN_OUTEDGE(head); (u=OUTVAL(e))!=0; e=NEXT_OUTEDGE(e)) { 
      if (u != tail){
        L2tu=ochange; /* L2tu will be # shrd prtnrs of (tail,u) not incl. head */
        /* step through inedges of u, incl. (head,u) itself */
        for(f = MIN_INEDGE(u); (v=INVAL(f))!=0; f=NEXT_INEDGE(f)) {
          if(IS_OUTEDGE(tail,v)) L2tu++;
        }
        cumchange += pow(oneexpa,(double)L2tu); /* sign corrected below */
      }
    }
    /* step through inedges of tail */
    for(e = MIN_INEDGE(tail); (u=INVAL(e))!=0; e=NEXT_INEDGE(e)) {
      if (u != head){
        L2uh=ochange; /* L2uh will be # shrd prtnrs of (u,head) not incl. tail */
        /* step through outedges of u , incl. (u,tail) itself */
        for(f = MIN_OUTEDGE(u);(v=OUTVAL(f))!=0; f=NEXT_OUTEDGE(f)){
          if(IS_OUTEDGE(v,head)) L2uh++;
        }
        cumchange += pow(oneexpa,(double)L2uh); /* sign corrected below */
      }
    }
    CHANGE_STAT[0] += echange * cumchange;
    TOGGLE_IF_MORE_TO_COME(i);
  }
  
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_gwtesp
*****************/
D_CHANGESTAT_FN(d_gwtesp) { 
  Edge e, f;
  int i, echange, ochange;
  int L2th, L2tu, L2uh;
  Vertex tail, head, u, v;
  double alpha, oneexpa, cumchange;
  
  CHANGE_STAT[0] = 0.0;
  alpha = INPUT_PARAM[0];
  oneexpa = 1.0-exp(-alpha);
  
  /* *** don't forget tail -> head */    
  FOR_EACH_TOGGLE(i){      
    cumchange=0.0;
    L2th=0;
    ochange = (EdgetreeSearch(tail=TAIL(i), head=HEAD(i), nwp->outedges) == 0) ? 0 : -1;
    echange = 2*ochange + 1;
    /* step through outedges of head  */
    for(e = EdgetreeMinimum(nwp->outedges, head);
	(u = nwp->outedges[e].value) != 0;
	e = EdgetreeSuccessor(nwp->outedges, e)){
      if (EdgetreeSearch(tail, u, nwp->outedges) != 0){
	L2tu=ochange;
	/* step through inedges of u */
	for(f = EdgetreeMinimum(nwp->inedges, u); 
	    (v = nwp->inedges[f].value) != 0;
	    f = EdgetreeSuccessor(nwp->inedges, f)){
	  if(EdgetreeSearch(tail,v,nwp->outedges)!= 0) L2tu++;
	}
	cumchange += pow(oneexpa,(double)L2tu);
      }
    }
    /* step through inedges of head */
    
    for(e = EdgetreeMinimum(nwp->inedges, head);
	(u = nwp->inedges[e].value) != 0;
	e = EdgetreeSuccessor(nwp->inedges, e)){
      if (EdgetreeSearch(tail, u, nwp->outedges) != 0){
	L2th++;
      }
      if (EdgetreeSearch(u, tail, nwp->outedges) != 0){
	L2uh=ochange;
	/* step through outedges of u */
	for(f = EdgetreeMinimum(nwp->outedges, u);
	    (v = nwp->outedges[f].value) != 0;
	    f = EdgetreeSuccessor(nwp->outedges, f)){
	  if(EdgetreeSearch(v,head,nwp->outedges)!= 0) L2uh++;
	}
	cumchange += pow(oneexpa,(double)L2uh) ;
      }
    }
    
    if(alpha < 100.0){
      cumchange += exp(alpha)*(1.0-pow(oneexpa,(double)L2th)) ;
    }else{
      cumchange += (double)L2th;
    }
    cumchange  = echange*cumchange;
    (CHANGE_STAT[0]) += cumchange;
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_gwtnsp
*****************/
D_CHANGESTAT_FN(d_gwtnsp) { 
  Edge e, f;
  int i, echange, ochange;
  int L2th, L2tu, L2uh;
  Vertex tail, head, u, v;
  double alpha, oneexpa, cumchange;
  
  CHANGE_STAT[0] = 0.0;

  alpha = INPUT_PARAM[0];
  oneexpa = 1.0-exp(-alpha);
  
  /* *** don't forget tail -> head */    
  FOR_EACH_TOGGLE(i){
    tail=TAIL(i); head=HEAD(i);
    cumchange=0.0;
    ochange = -IS_OUTEDGE(tail,head);
    echange = 2*ochange + 1;
    /* step through outedges of head */
    for(e = MIN_OUTEDGE(head); (u=OUTVAL(e))!=0; e=NEXT_OUTEDGE(e)) { 
      if (u != tail){
        L2tu=ochange; /* L2tu will be # shrd prtnrs of (tail,u) not incl. head */
        /* step through inedges of u, incl. (head,u) itself */
        for(f = MIN_INEDGE(u); (v=INVAL(f))!=0; f=NEXT_INEDGE(f)) {
          if(IS_OUTEDGE(tail,v)) L2tu++;
        }
        cumchange += pow(oneexpa,(double)L2tu); /* sign corrected below */
      }
    }
    /* step through inedges of tail */
    for(e = MIN_INEDGE(tail); (u=INVAL(e))!=0; e=NEXT_INEDGE(e)) {
      if (u != head){
        L2uh=ochange; /* L2uh will be # shrd prtnrs of (u,head) not incl. tail */
        /* step through outedges of u , incl. (u,tail) itself */
        for(f = MIN_OUTEDGE(u);(v=OUTVAL(f))!=0; f=NEXT_OUTEDGE(f)){
          if(IS_OUTEDGE(v,head)) L2uh++;
        }
        cumchange += pow(oneexpa,(double)L2uh); /* sign corrected below */
      }
    }
    CHANGE_STAT[0] += echange * cumchange;
    TOGGLE_IF_MORE_TO_COME(i);
  }
  
  UNDO_PREVIOUS_TOGGLES(i);

  alpha = INPUT_PARAM[0];
  oneexpa = 1.0-exp(-alpha);
  
  /* *** don't forget tail -> head */    
  FOR_EACH_TOGGLE(i){      
    cumchange=0.0;
    L2th=0;
    ochange = (EdgetreeSearch(tail=TAIL(i), head=HEAD(i), nwp->outedges) == 0) ? 0 : -1;
    echange = 2*ochange + 1;
    /* step through outedges of head  */
    for(e = EdgetreeMinimum(nwp->outedges, head);
	(u = nwp->outedges[e].value) != 0;
	e = EdgetreeSuccessor(nwp->outedges, e)){
      if (EdgetreeSearch(tail, u, nwp->outedges) != 0){
	L2tu=ochange;
	/* step through inedges of u */
	for(f = EdgetreeMinimum(nwp->inedges, u); 
	    (v = nwp->inedges[f].value) != 0;
	    f = EdgetreeSuccessor(nwp->inedges, f)){
	  if(EdgetreeSearch(tail,v,nwp->outedges)!= 0) L2tu++;
	}
	cumchange += pow(oneexpa,(double)L2tu);
      }
    }
    /* step through inedges of head */
    
    for(e = EdgetreeMinimum(nwp->inedges, head);
	(u = nwp->inedges[e].value) != 0;
	e = EdgetreeSuccessor(nwp->inedges, e)){
      if (EdgetreeSearch(tail, u, nwp->outedges) != 0){
	L2th++;
      }
      if (EdgetreeSearch(u, tail, nwp->outedges) != 0){
	L2uh=ochange;
	/* step through outedges of u */
	for(f = EdgetreeMinimum(nwp->outedges, u);
	    (v = nwp->outedges[f].value) != 0;
	    f = EdgetreeSuccessor(nwp->outedges, f)){
	  if(EdgetreeSearch(v,head,nwp->outedges)!= 0) L2uh++;
	}
	cumchange += pow(oneexpa,(double)L2uh) ;
      }
    }
    
    if(alpha < 100.0){
      cumchange += exp(alpha)*(1.0-pow(oneexpa,(double)L2th)) ;
    }else{
      cumchange += (double)L2th;
    }
    cumchange  = echange*cumchange;
    (CHANGE_STAT[0]) -= cumchange;
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);

}

/********************  changestats:  H    ***********/
/*****************
 changestat: d_hamming
*****************/
/* This function must be passed two networks in the forms of edgelists:
   One is the network from which hamming distances are calculated and the
   other is the network of weights on the dyads, which is an edgelist with a
   third column of weights.  Note that all non-edges in the second network
   will have the value given by defaultval; thus, the "unweighted" hamming
   distance is obtained when the default is 1.0 and the second network is
   empty. */
D_CHANGESTAT_FN(d_hamming) { 
  int i, discord;
  
  ZERO_ALL_CHANGESTATS(i);
  Edge wt_net_start= INPUT_PARAM[0]*2+2;
  double defaultval = INPUT_PARAM[wt_net_start-1]; /* Hamming wt for non-edges in cov nw */
  double *wt_net = INPUT_PARAM+wt_net_start;

  /* *** don't forget tail -> head */    
  FOR_EACH_TOGGLE(i) {
    Vertex tail=TAIL(i), head=HEAD(i);
    
    discord = XOR(dEdgeListSearch(tail, head, INPUT_PARAM), IS_OUTEDGE(tail, head));

    /* Second, search second network to see if the weight is different from
       defaultval.  In unweighted case, this network is empty. */

    Edge wt_pos = dEdgeListSearch(tail, head, wt_net);
    double val = wt_pos ? wt_net[wt_pos+2*(unsigned int)wt_net[0]] : defaultval;

    CHANGE_STAT[0] += (discord ? -val : val);
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);  
}

/*****************
 changestat: d_hammingmix_constant
*****************/
D_CHANGESTAT_FN(d_hammingmix_constant) { 
  int i, nhedge, discord;
  int matchvaltail, matchvalhead;
  
  nhedge = INPUT_PARAM[0];
/*  Rprintf("nhedge %d\n", nhedge); */
  if(ntoggles==2){
   matchvaltail = INPUT_PARAM[TAIL(0)+2*nhedge];
   matchvalhead = INPUT_PARAM[HEAD(0)+2*nhedge];
   if(matchvaltail != INPUT_PARAM[TAIL(1)+2*nhedge] ||
      matchvalhead != INPUT_PARAM[HEAD(1)+2*nhedge]){
      CHANGE_STAT[0] = 10000.0;
      return;
   }
  }
  CHANGE_STAT[0] = 0.0;
/*  Rprintf("Warning: hammingconstantmix can only be used with ConstantEdges terms.\n");
  Rprintf("nhedge %d i0 %f i1 %f i2 %f i3 %f\n", nhedge, INPUT_PARAM[0],
                                 INPUT_PARAM[1],
                                 INPUT_PARAM[2],
                                 INPUT_PARAM[3]
		  ); */
     

  /* *** don't forget tail -> head */    
  FOR_EACH_TOGGLE(i)
    {
      discord=XOR(dEdgeListSearch(TAIL(i), HEAD(i), INPUT_PARAM), IS_OUTEDGE(TAIL(i), HEAD(i)));
      CHANGE_STAT[0] += (discord ? -1.0 : 1.0);

    if (i+1 < ntoggles){
      ToggleEdge(TAIL(i), HEAD(i), &nwp[0]);  /* Toggle this edge if more to come */
    }
  }
  i--;
  while (--i>=0){  /*  Undo all previous toggles. */
    ToggleEdge(TAIL(i), HEAD(i), &nwp[0]);
  }
}

/*****************
 changestat: d_hammingmix
*****************/
D_CHANGESTAT_FN(d_hammingmix) { 
  int i;
  
  Edge nhedge =  INPUT_PARAM[0];
/*  Rprintf("nstats %d nhedge %d i0 %f i1 %f i2 %f i3 %f\n",nstats, nhedge, INPUT_PARAM[0],
                                 INPUT_PARAM[1],
                                 INPUT_PARAM[2],
                                 INPUT_PARAM[3]
		  ); */


  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i){
    Vertex tail=TAIL(i), head=HEAD(i);
    int matchvaltail = INPUT_PARAM[tail+2*N_CHANGE_STATS+2*nhedge];
    int matchvalhead = INPUT_PARAM[head+2*N_CHANGE_STATS+2*nhedge];
    unsigned int discord = XOR(dEdgeListSearch(tail, head, INPUT_PARAM), IS_OUTEDGE(tail, head));
    for (unsigned int j=0; j<N_CHANGE_STATS; j++){
      if(matchvaltail==INPUT_PARAM[2*nhedge+1+j] &&
        matchvalhead==INPUT_PARAM[2*nhedge+1+N_CHANGE_STATS+j]){
          CHANGE_STAT[j] += (discord ? -1.0 : 1.0);
      }
    }
      
    if (i+1 < ntoggles){
      ToggleEdge(TAIL(i), HEAD(i), &nwp[0]);  /* Toggle this edge if more to come */
    }
  }
  i--;
  while (--i>=0){  /*  Undo all previous toggles. */
    ToggleEdge(TAIL(i), HEAD(i), &nwp[0]);
  }
}

/********************  changestats:  I    ***********/

// A macro indicating whether x is in [from,to)
#define FROM_TO(x, from, to) ((x)>=(from) && (x)<(to))

/*****************
 changestat: d_idegrange
*****************/
D_CHANGESTAT_FN(d_idegrange) { 
  int i, j, echange;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    Vertex tail, head;
    echange=(EdgetreeSearch(tail=TAIL(i), head=HEAD(i), nwp->outedges)==0)? 1:-1;
    Vertex headideg = IN_DEG[head];
    for(j = 0; j < N_CHANGE_STATS; j++) {
      Vertex from = INPUT_PARAM[2*j], to = INPUT_PARAM[2*j+1];
      CHANGE_STAT[j] += FROM_TO(headideg + echange, from, to) - FROM_TO(headideg, from, to);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}
 
/*****************
 changestat: d_idegrange_by_attr
*****************/
D_CHANGESTAT_FN(d_idegrange_by_attr) { 
  /* The inputparams are assumed to be set up as follows:
  The first 3*nstats values are in triples:  (from, to, attrvalue)
  The values following the first 2*nstats values are the nodal attributes.
  */
  int i, j;
  Vertex *id;
  
  id=IN_DEG;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    Vertex tail, head;
    int echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1:1;
    Vertex headideg = id[head];
    int headattr = INPUT_PARAM[3*N_CHANGE_STATS + head - 1]; 
    for(j = 0; j < N_CHANGE_STATS; j++){
      Vertex from = INPUT_PARAM[3*j], to = INPUT_PARAM[3*j + 1];
      int testattr = INPUT_PARAM[3*j + 2]; 
      if (headattr == testattr)  /* we have tail attr match */
        CHANGE_STAT[j] += FROM_TO(headideg + echange, from, to) - FROM_TO(headideg, from, to);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_idegrange_w_homophily
*****************/
D_CHANGESTAT_FN(d_idegrange_w_homophily) { 
  /*  The inputparams are assumed to be set up as follows:
  The first 2*nstats values are the values of idegrange
  The values following the first 2*nstats values are the nodal attributes.
  */
  int i, j;
  double *nodeattr;
  Edge e;

  nodeattr = mtp->inputparams + N_CHANGE_STATS*2 - 1;  

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    Vertex tail=TAIL(i), head=HEAD(i);
    int tailattr = nodeattr[tail], headattr = nodeattr[head];
    if (headattr == tailattr) { /* They match; otherwise don't bother */
      int echange = IS_OUTEDGE(tail, head) ? -1:1;
      Vertex headideg=0, v;
      STEP_THROUGH_INEDGES(head, e, v) { headideg += (nodeattr[v]==headattr); }
      for(j = 0; j < N_CHANGE_STATS; j++) {
        Vertex from = INPUT_PARAM[2*j], to = INPUT_PARAM[2*j+1];
        CHANGE_STAT[j] += FROM_TO(headideg + echange, from, to) - FROM_TO(headideg, from, to);
      }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}                                        

#undef FROM_TO

/*****************
 changestat: d_idegree
*****************/
D_CHANGESTAT_FN(d_idegree) { 
  int i, j;
  
  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i){
    Vertex head;
    int echange = (EdgetreeSearch(TAIL(i), head=HEAD(i), nwp->outedges) == 0) ? 1 : -1;
    Vertex headd = IN_DEG[head];
    
    for(j=0; j < N_CHANGE_STATS; j++){
      Vertex deg = INPUT_PARAM[j];
      CHANGE_STAT[j] += (headd + echange == deg) - (headd == deg);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_idegree_by_attr
*****************/
D_CHANGESTAT_FN(d_idegree_by_attr) { 
  /* The inputparams are assumed to be set up as follows:
  The first 2*nstats values are in pairs:  (degree, attrvalue)
  The values following the first 2*nstats values are the nodal attributes.
  */
  int i, j, echange, headattr, testattr;
  Vertex head, headdeg, d, *id;
  
  id=IN_DEG;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    echange=(EdgetreeSearch(TAIL(i), head=HEAD(i), nwp->outedges)==0)? 1:-1;
    headdeg = id[head];
    headattr = INPUT_PARAM[2*N_CHANGE_STATS + head - 1]; 
    for(j = 0; j < N_CHANGE_STATS; j++) {
      d = (Vertex)INPUT_PARAM[2*j];
      testattr = INPUT_PARAM[2*j + 1]; 
      if (headattr == testattr)  /* we have head attr match */
        CHANGE_STAT[j] += (headdeg + echange == d) - (headdeg == d);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_idegree_w_homophily
*****************/
D_CHANGESTAT_FN(d_idegree_w_homophily) { 
  /*  The inputparams are assumed to be set up as follows:
  The first nstats values are the values of degree
  The values following the first nstats values are the nodal attributes.
  */
  int i, j, echange, tailattr, headattr;
  Vertex tail, head, headdeg, deg, tmp;
  double *nodeattr;
  Edge e;

  nodeattr = mtp->inputparams + N_CHANGE_STATS - 1;  

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    tail=TAIL(i);
    head=HEAD(i);
    tailattr = (int)nodeattr[tail];
    headattr = (int)nodeattr[head];    
    if (tailattr == headattr) { /* They match; otherwise don't bother */
      echange=(EdgetreeSearch(tail, head, nwp->outedges)==0)? 1:-1;
      headdeg=0;
/*      for(e = EdgetreeMinimum(nwp->outedges, head);
      (tmp = nwp->outedges[e].value) != 0;
      e = EdgetreeSuccessor(nwp->outedges, e)) {
        headdeg += (nodeattr[tmp]==headattr);
      } */
      STEP_THROUGH_INEDGES(head, e, tmp){
        headdeg += (nodeattr[tmp]==headattr);
      }
      for(j = 0; j < N_CHANGE_STATS; j++) {
        deg = (Vertex)INPUT_PARAM[j];
        CHANGE_STAT[j] += (headdeg + echange == deg) - (headdeg == deg);
      }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_idegreepopularity
*****************/
D_CHANGESTAT_FN(d_idegreepopularity) { 
  int i, edgeflag;
  double change;
  Vertex head, tail, deg=0;
  
  /* *** don't forget tail -> head */    
  change = 0.0;
  FOR_EACH_TOGGLE(i) {
    tail=TAIL(i);
    head=HEAD(i);
    edgeflag = IS_OUTEDGE(tail, head); /* either 0 or 1 */
    deg = (double)(IN_DEG[head]);
    if(edgeflag){
      change -= sqrt(deg);
      change += (deg-1.0)*(sqrt(deg-1.0)-sqrt(deg));
    }else{
      change += sqrt(deg+1.0);
      change += deg*(sqrt(deg+1.0)-sqrt(deg));
    }
    TOGGLE_IF_MORE_TO_COME(i); 
  }
  CHANGE_STAT[0]=change; 
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_intransitive
*****************/
D_CHANGESTAT_FN(d_intransitive) { 
  Edge e;
  Vertex tail, head, node2;
  double change;
  int edgeflag, i;
  
  /* *** don't forget tail -> head */    
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i) 
  {
    edgeflag = IS_OUTEDGE(tail = TAIL(i), head = HEAD(i));
    change = 0.0;
    STEP_THROUGH_OUTEDGES(head, e, node2) {
      if (node2 != tail){
        if (EdgetreeSearch(tail, node2, nwp->outedges) == 0){
          change = change + 1.0;
        }
      }
    }
    STEP_THROUGH_INEDGES(head, e, node2) {
      if (node2 != tail){
        if (EdgetreeSearch(tail, node2, nwp->outedges) != 0){
          change = change - 1.0;
        }
      }
    }
    STEP_THROUGH_INEDGES(tail, e, node2) {
      if (node2 != head){
        if (EdgetreeSearch(node2, head, nwp->outedges) == 0){
          change = change + 1.0;
        }
      }
    }    
    CHANGE_STAT[0] += edgeflag ? -change : change;
/*  Rprintf("tail %d head %d edgeflag %d change %f\n",tail,head, edgeflag, change); */
    TOGGLE_IF_MORE_TO_COME(i);
  }  
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_isolates
*****************/
D_CHANGESTAT_FN(d_isolates) { 
  int i, echange;
  Vertex tail, head, taild, headd=0, *id, *od;

  id=IN_DEG;
  od=OUT_DEG;
  
  /* *** don't forget tail -> head */    
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i)
    {      
      echange = (EdgetreeSearch(tail=TAIL(i), head=HEAD(i), nwp->outedges) == 0) ? 1 : -1;
      taild = od[tail] + id[tail];
      headd = od[head] + id[head];
      CHANGE_STAT[0] += (taild + echange == 0) - (taild == 0);
      CHANGE_STAT[0] += (headd + echange == 0) - (headd == 0);
      
      TOGGLE_IF_MORE_TO_COME(i);
    }

  UNDO_PREVIOUS_TOGGLES(i);
}

S_CHANGESTAT_FN(s_isolates) { 
  /* int i, echange;
     Vertex tail, head, taild, headd=0, *id, *od; */
  Vertex *id, *od;

  id=IN_DEG;
  od=OUT_DEG;
  
  /* *** don't forget tail -> head */    
  CHANGE_STAT[0] = 0.0;
  for(Vertex tail=1; tail <= N_NODES; tail++){
    if(od[tail] + id[tail] == 0)
      CHANGE_STAT[0] ++;
  }
}

/*****************
 changestat: d_istar
*****************/
D_CHANGESTAT_FN(d_istar) { 
  double change, headd=0.0;
  int edgeflag, i, j, kmo;
  Edge e;
  Vertex tail, head, node3;
  int ninputs, nstats;
  double tailattr;
  
  ninputs = (int)N_INPUT_PARAMS;
  nstats  = (int)N_CHANGE_STATS;
  
  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  if(ninputs>nstats){
    /* match on attributes */
    for (i=0; i < ntoggles; i++) {
      /* edgeflag is 1 if edge exists and will disappear
      edgeflag is 0 if edge DNE and will appear */
      edgeflag = IS_OUTEDGE(tail = TAIL(i), head = HEAD(i));
      tailattr = INPUT_ATTRIB[tail-1];
      if(tailattr == INPUT_ATTRIB[head-1]){
        headd = - edgeflag;
        for(e = EdgetreeMinimum(nwp->inedges, head);
	      (node3 = nwp->inedges[e].value) != 0;
	      e = EdgetreeSuccessor(nwp->inedges, e)) {/* step through inedges of head */
          if(tailattr == INPUT_ATTRIB[node3-1]){++headd;}
        }	  
        for(j=0; j < N_CHANGE_STATS; j++) {
          kmo = ((int)INPUT_PARAM[j]) - 1;
          change = CHOOSE(headd, kmo); 
          CHANGE_STAT[j] += (edgeflag ? - change : change); 
        }
      }
      TOGGLE_IF_MORE_TO_COME(i);
    }
  }else{
    for (i=0; i < ntoggles; i++) {
      /* edgeflag is 1 if edge exists and will disappear
      edgeflag is 0 if edge DNE and will appear */
      edgeflag = IS_OUTEDGE(tail = TAIL(i), head = HEAD(i));
      headd = IN_DEG[head] - edgeflag;	
      for(j=0; j < N_CHANGE_STATS; j++) {
        kmo = ((int)INPUT_PARAM[j]) - 1;
        change = CHOOSE(headd, kmo); 
        CHANGE_STAT[j] += (edgeflag ? - change : change); 
      }
      TOGGLE_IF_MORE_TO_COME(i);
    }
  }
  
  UNDO_PREVIOUS_TOGGLES(i);
}

/********************  changestats:  K    ***********/
/*****************
 changestat: d_kstar
*****************/
D_CHANGESTAT_FN(d_kstar) { 
  double change, taild, headd=0.0;
  int edgeflag, i, j, kmo;
  Edge e;
  Vertex tail, head, node3;
  int ninputs, nstats;
  double tailattr;
    
  ninputs = (int)N_INPUT_PARAMS;
  nstats  = (int)N_CHANGE_STATS;
  
  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  if(ninputs>nstats){
    /* match on attributes */
    for (i=0; i < ntoggles; i++) {
      /* edgeflag is 1 if edge exists and will disappear
      edgeflag is 0 if edge DNE and will appear */
      edgeflag = IS_OUTEDGE(tail = TAIL(i), head = HEAD(i));
      tailattr = INPUT_ATTRIB[tail-1];
      if(tailattr == INPUT_ATTRIB[head-1]){
        taild = - edgeflag;
        STEP_THROUGH_OUTEDGES(tail, e, node3) {
          if(tailattr == INPUT_ATTRIB[node3-1]){++taild;}
        }
        STEP_THROUGH_INEDGES(tail, e, node3) {
          if(tailattr == INPUT_ATTRIB[node3-1]){++taild;}
        }
        headd = - edgeflag;
        STEP_THROUGH_OUTEDGES(head, e, node3) {
          if(tailattr == INPUT_ATTRIB[node3-1]){++headd;}
        }
        STEP_THROUGH_INEDGES(head, e, node3) {
          if(tailattr == INPUT_ATTRIB[node3-1]){++headd;}
        }
        
        for(j=0; j < N_CHANGE_STATS; j++) {
          kmo = ((int)INPUT_PARAM[j]) - 1;
/*          if (kmo==0) {
            change=1;
          } else { */
            change = CHOOSE(taild, kmo) + CHOOSE(headd, kmo); 
/*          } uncomment these few lines to define 1-stars as equivalent to 
              edges (currently, each edge is counted as two 1-stars) */
          CHANGE_STAT[j] += (edgeflag ? - change : change); 
        }
      }
      TOGGLE_IF_MORE_TO_COME(i);
    }
  }else{
    /* *** don't forget tail -> head */    
    for (i=0; i < ntoggles; i++)
    {
      /* edgeflag is 1 if edge exists and will disappear
      edgeflag is 0 if edge DNE and will appear */
      edgeflag = IS_OUTEDGE(tail = TAIL(i), head = HEAD(i));
      taild = OUT_DEG[tail] + IN_DEG[tail] - edgeflag; 
      headd = OUT_DEG[head] + IN_DEG[head] - edgeflag;
      for(j=0; j < N_CHANGE_STATS; j++) 
      {
        kmo = ((int)INPUT_PARAM[j]) - 1;
/*        if (kmo==0) {
          change=1;
        } else { */
          change = CHOOSE(taild, kmo) + CHOOSE(headd, kmo); 
/*      } uncomment these few lines to define 1-stars as equivalent to 
          edges (currently, each edge is counted as two 1-stars) */
        CHANGE_STAT[j] += (edgeflag ? - change : change); 
      }
      TOGGLE_IF_MORE_TO_COME(i);
    }
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/********************  changestats:  L    ***********/
/*****************
 changestat: d_localtriangle
*****************/
D_CHANGESTAT_FN(d_localtriangle) { 
  Edge e;
  Vertex tail, head, node3, nmat;
  double change;
  int edgeflag, i;
  
  nmat = (Vertex)(INPUT_PARAM[0]);
  
  /* *** don't forget tail -> head */    
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i) 
    {
      edgeflag = IS_OUTEDGE(tail = TAIL(i), head = HEAD(i));
      change = 0.0;
      
      if(INPUT_PARAM[1+(HEAD(i)-1)+(TAIL(i)-1)*nmat] == 1.0){
        STEP_THROUGH_OUTEDGES(head, e, node3) { /* step through outedges of head */
	    if(INPUT_PARAM[1+(node3-1)+(TAIL(i)-1)*nmat] == 1.0 && 
	       INPUT_PARAM[1+(node3-1)+(HEAD(i)-1)*nmat] == 1.0 ){
	      if (DIRECTED){
		if (IS_INEDGE(node3,tail) ) ++change;
		if (IS_OUTEDGE(node3,tail)) ++change;
	      }else{
		if (IS_UNDIRECTED_EDGE(node3,tail)) ++change;
	      }
	    }
	  }
	
        STEP_THROUGH_INEDGES(head, e, node3) { /* step through inedges of head */
	    if(INPUT_PARAM[1+(node3-1)+(TAIL(i)-1)*nmat] == 1.0 && 
	       INPUT_PARAM[1+(node3-1)+(HEAD(i)-1)*nmat] == 1.0 ){
	      if (DIRECTED)
		{
		if (IS_INEDGE(node3,tail) ) ++change;
		if (IS_OUTEDGE(node3,tail)) ++change;
		}
	      else
		{
		  if (IS_UNDIRECTED_EDGE(node3,tail)) ++change;
		}
	    }
	  }
	
	CHANGE_STAT[0] += edgeflag ? - change : change;
      
      }
      TOGGLE_IF_MORE_TO_COME(i);
    }  
  UNDO_PREVIOUS_TOGGLES(i);
}

/********************  changestats:  M    ***********/
/*****************
 changestat: d_m2star
*****************/
D_CHANGESTAT_FN(d_m2star) {
  Vertex tail, head;
  int tailid, headod, change;
  int i, edgeflag, backedgeflag;
    
  CHANGE_STAT[0] = 0.0;

  /* *** don't forget tail -> head */    
  for (i=0; i < ntoggles; i++)
    {
      /*  edgeflag is 1 if the edge from TAIL(i) to HEAD(i)  */
      /*   exists and will disappear */
      /*  edgeflag is 0 if the edge does not exist */
      edgeflag = IS_OUTEDGE(tail = TAIL(i), head = HEAD(i));
      backedgeflag = (EdgetreeSearch(head, tail, nwp->outedges) != 0);

      tailid = IN_DEG[tail]; 
      headod = OUT_DEG[head];
      change = tailid + headod - 2*backedgeflag; 
      CHANGE_STAT[0] += (edgeflag ? -change : change); 

      TOGGLE_IF_MORE_TO_COME(i);
    }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_meandeg
*****************/
D_CHANGESTAT_FN(d_meandeg) {
  int i;

  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i) {
    CHANGE_STAT[0] += (IS_OUTEDGE(TAIL(i), HEAD(i)) ? -2.0 : 2.0);
    TOGGLE_IF_MORE_TO_COME(i);
  }
  CHANGE_STAT[0]/=(double)N_NODES;
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_mix
 This appears to be the version of nodemix used for 
 bipartite networks (only)
*****************/
D_CHANGESTAT_FN(d_mix) {
  Vertex tail, head, tmpi;
  int matchvaltail, matchvalhead;
  int i, j, edgeflag, nstats;

  nstats = N_CHANGE_STATS;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    tail=TAIL(i);
    head=HEAD(i);
    edgeflag = IS_OUTEDGE(tail, head);
    if (BIPARTITE > 0 && tail > head) { 
      tmpi = tail; tail = head; head = tmpi; /* swap tail, head */
    }
    matchvaltail = INPUT_PARAM[tail-1+2*nstats];
    matchvalhead = INPUT_PARAM[head-1+2*nstats];
    for (j=0; j<nstats; j++) {
      if(matchvaltail==INPUT_PARAM[j] && matchvalhead==INPUT_PARAM[nstats+j]) {
        CHANGE_STAT[j] += edgeflag ? -1.0 : 1.0;
      }
	  }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_mutual

 (1,1) -> anything = -1
 anything -> (1,1) = +1
*****************/
D_CHANGESTAT_FN(d_mutual) { 
  double matchval, change;
  Vertex tail, head;
  int i, j, ninputs, noattr;

  ninputs = N_INPUT_PARAMS - N_NODES;
  noattr = (N_INPUT_PARAMS == 0);

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    tail = TAIL(i);
    head = HEAD(i);
    if (IS_OUTEDGE(head,tail)) { /* otherwise, no change occurs */
      change = IS_OUTEDGE(tail, head) ? -1.0 : 1.0 ;
      if (noattr) { /* "plain vanilla" mutual, without node attributes */
        CHANGE_STAT[0] += change;
      } else { /* Only consider mutuals where node attributes match */
        matchval = INPUT_PARAM[tail+ninputs-1];
        if (matchval == INPUT_PARAM[head+ninputs-1]) { /* We have a match! */
          if (ninputs==0) {/* diff=F in network statistic specification */
            CHANGE_STAT[0] += change;
          } else { /* diff=T */
            for (j=0; j<ninputs; j++) {
              if (matchval == INPUT_PARAM[j]) 
                CHANGE_STAT[j] += change;
            }
          }
        }
      }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_mutual_by_attr
*****************/
D_CHANGESTAT_FN(d_mutual_by_attr) { 
  double change;
  Vertex tail, head;
  int i, j, ninputs;

  ninputs = N_INPUT_PARAMS - N_NODES;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    tail = TAIL(i);
    head = HEAD(i);
    if (IS_OUTEDGE(head,tail)) { /* otherwise, no change occurs */
      change = IS_OUTEDGE(tail, head) ? -1.0 : 1.0 ;
      for (j=0; j<ninputs; j++) {
        if (INPUT_PARAM[tail+ninputs-1] == INPUT_PARAM[j]){CHANGE_STAT[j] += change;}
        if (INPUT_PARAM[head+ninputs-1] == INPUT_PARAM[j]){CHANGE_STAT[j] += change;}
      }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/********************  changestats:  N    ***********/
/*****************
 changestat: d_nearsimmelian
*****************/
D_CHANGESTAT_FN(d_nearsimmelian) { 
  Vertex tail, head, node3;
  double change;
  int edgeflag, i, edgeflagth, sc;

  /* *** don't forget tail -> head */    
 CHANGE_STAT[0] = 0.0;
 FOR_EACH_TOGGLE(i) {
  edgeflag = IS_OUTEDGE(tail = TAIL(i), head = HEAD(i));
  edgeflagth = (EdgetreeSearch(head, tail, nwp->outedges) == 0);
   
  for(node3=1;node3<=N_NODES;node3++){
    if((node3!=tail)&&(node3!=head)){
     sc = edgeflagth + (EdgetreeSearch(node3, tail, nwp->outedges) == 0);
     if(sc < 2){
      sc += (EdgetreeSearch(tail, node3, nwp->outedges) == 0);
      if(sc < 2){
       sc += (EdgetreeSearch(node3, head, nwp->outedges) == 0);
       if(sc < 2){
        sc += (EdgetreeSearch(head, node3, nwp->outedges) == 0);
        if(sc < 2){
         change=0.0;
         if (sc == 0 && edgeflag == 0 ){--change;}
         if (sc == 0 && edgeflag == 1 ){++change;}
         if (sc == 1 && edgeflag == 0 ){++change;}
         if (sc == 1 && edgeflag == 1 ){--change;}
         CHANGE_STAT[0] += change;
	}
       }
      }
     }
    }
   }
   
   TOGGLE_IF_MORE_TO_COME(i);
  }
  
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_nodecov
*****************/
D_CHANGESTAT_FN(d_nodecov) { 
  double sum;
  Vertex tail, head;
  int i, edgeflag;

  /* *** don't forget tail -> head */    
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i) 
    {
      edgeflag=IS_OUTEDGE(tail = TAIL(i), head = HEAD(i));
      sum = INPUT_ATTRIB[tail-1] + INPUT_ATTRIB[head-1];
      CHANGE_STAT[0] += edgeflag ? -sum : sum;
      TOGGLE_IF_MORE_TO_COME(i);
    }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_nodefactor
*****************/
D_CHANGESTAT_FN(d_nodefactor) { 
  double s, factorval;
  Vertex tail, head;
  int i, j, tailattr, headattr;
  
  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    tail = TAIL(i);
    head = HEAD(i);
    s = IS_OUTEDGE(tail, head) ? -1.0 : 1.0;
    tailattr = INPUT_ATTRIB[tail-1];
    headattr = INPUT_ATTRIB[head-1];
    for (j=0; j < N_CHANGE_STATS; j++) {
      factorval = INPUT_PARAM[j];
      if (tailattr == factorval) CHANGE_STAT[j] += s;
      if (headattr == factorval) CHANGE_STAT[j] += s;
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_nodeicov
*****************/
D_CHANGESTAT_FN(d_nodeicov) { 
  double sum;
  Vertex tail, head;
  int i, edgeflag;
  
  /* *** don't forget tail -> head */    
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i) 
    {
      edgeflag=IS_OUTEDGE(tail = TAIL(i), head = HEAD(i));
      sum = INPUT_ATTRIB[head-1];
      CHANGE_STAT[0] += edgeflag ? -sum : sum;
      TOGGLE_IF_MORE_TO_COME(i);
    }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_nodeifactor
*****************/
D_CHANGESTAT_FN(d_nodeifactor) { 
  double s;
  Vertex head;
  int i, j, headattr;
  
  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    head = HEAD(i);
    s = IS_OUTEDGE(TAIL(i), head) ? -1.0 : 1.0;
    headattr = INPUT_ATTRIB[head-1];
    for (j=0; j < N_CHANGE_STATS; j++) {
      if (headattr == INPUT_PARAM[j]) CHANGE_STAT[j] += s;
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_nodematch
*****************/
D_CHANGESTAT_FN(d_nodematch) { 
  double matchval;
  Vertex tail, head, ninputs;
  int i, j, edgeflag;
  
  ninputs = N_INPUT_PARAMS - N_NODES;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    tail=TAIL(i);
    head=HEAD(i);
    matchval = INPUT_PARAM[tail+ninputs-1];
    if (matchval == INPUT_PARAM[head+ninputs-1]) { /* We have a match! */
      edgeflag = IS_OUTEDGE(tail, head);
      if (ninputs==0) {/* diff=F in network statistic specification */
        CHANGE_STAT[0] += edgeflag ? -1.0 : 1.0;
      } else { /* diff=T */
        for (j=0; j<ninputs; j++) {
          if (matchval == INPUT_PARAM[j]) 
            CHANGE_STAT[j] += edgeflag ? -1.0 : 1.0;
        }
      }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_nodemix
 Update mixing matrix, non-bipartite networks only 
 (but see also d_mix)
*****************/
D_CHANGESTAT_FN(d_nodemix) {
  Vertex tail, head;
  int i, j, ninputs, ninputs2;
  double rtype, ctype, tmp, change;

  ninputs = N_INPUT_PARAMS - N_NODES;
  ninputs2 = ninputs/2;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
      tail=TAIL(i);
      head=HEAD(i);
      change = IS_OUTEDGE(tail, head) ? -1.0 : 1.0;
      /*Find the node covariate values (types) for the tail and head*/
      rtype=INPUT_PARAM[tail+ninputs-1];
      ctype=INPUT_PARAM[head+ninputs-1];
      if (!DIRECTED && rtype > ctype)  {
        tmp = rtype; rtype = ctype; ctype = tmp; /* swap rtype, ctype */
      }
      /*Find the right statistic to update */
      for(j=0; j<ninputs2; j++){
        if((INPUT_PARAM[j] == rtype) && (INPUT_PARAM[j+ninputs2] == ctype)){
          CHANGE_STAT[j] += change;
          j = ninputs2; /* leave the for loop */
        }
      } 
      TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_nodeocov
*****************/
D_CHANGESTAT_FN(d_nodeocov) { 
  double sum;
  Vertex tail, head;
  int i, edgeflag;
  
  /* *** don't forget tail -> head */    
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i) 
    {
      edgeflag=IS_OUTEDGE(tail = TAIL(i), head = HEAD(i));
      sum = INPUT_ATTRIB[tail-1];
      CHANGE_STAT[0] += edgeflag ? -sum : sum;
      TOGGLE_IF_MORE_TO_COME(i);
    }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_nodeofactor
*****************/
D_CHANGESTAT_FN(d_nodeofactor) { 
  double s;
  Vertex tail;
  int i, j, tailattr;
  
  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    tail = TAIL(i);
    s = IS_OUTEDGE(tail, HEAD(i)) ? -1.0 : 1.0;
    tailattr = INPUT_ATTRIB[tail-1];
    for (j=0; j < N_CHANGE_STATS; j++) {
      if (tailattr == INPUT_PARAM[j]) CHANGE_STAT[j] += s;
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_nsp
*****************/
D_CHANGESTAT_FN(d_nsp) { 
  Edge e, f;
  int i, j, echange;
  int L2th, L2tu, L2uh;
  Vertex deg;
  Vertex tail, head, u, v;
  
  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 1;
    STEP_THROUGH_OUTEDGES(head, e, u) {
      if (u != tail){
        L2tu=0;
        STEP_THROUGH_OUTEDGES(u, f, v) {
          if(IS_OUTEDGE(MIN(v,tail),MAX(v,tail))) L2tu++;
        }
        STEP_THROUGH_INEDGES(u, f, v) {
          if(IS_OUTEDGE(MIN(v,tail),MAX(v,tail))) L2tu++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)INPUT_PARAM[j];
          CHANGE_STAT[j] += ((L2tu + echange == deg)
          - (L2tu == deg));
        }
      }
    }
    STEP_THROUGH_INEDGES(head, e, u) {
      if (u != tail){
        L2tu=0;
        STEP_THROUGH_OUTEDGES(u, f, v) {
          if(IS_OUTEDGE(MIN(v,tail),MAX(v,tail))) L2tu++;
        }
        STEP_THROUGH_INEDGES(u, f, v) {
          if(IS_OUTEDGE(MIN(v,tail),MAX(v,tail))) L2tu++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)INPUT_PARAM[j];
          CHANGE_STAT[j] += ((L2tu + echange == deg)
          - (L2tu == deg));
        }
      }
    }
    STEP_THROUGH_OUTEDGES(tail, e, u) {
      if (u != head){
        L2uh=0;
        STEP_THROUGH_OUTEDGES(u, f, v) {
          if(IS_OUTEDGE(MIN(v,head),MAX(v,head))) L2uh++;
        }
        STEP_THROUGH_INEDGES(u, f, v) {
          if(IS_OUTEDGE(MIN(v,head),MAX(v,head))) L2uh++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)INPUT_PARAM[j];
          CHANGE_STAT[j] += ((L2uh + echange == deg)
          - (L2uh == deg));
        }
      }
    }
    STEP_THROUGH_INEDGES(tail, e, u) {
      if (u != head){
        L2uh=0;
        STEP_THROUGH_OUTEDGES(u, f, v) {
          if(IS_OUTEDGE(MIN(v,head),MAX(v,head))) L2uh++;
        }
        STEP_THROUGH_INEDGES(u, f, v) {
          if(IS_OUTEDGE(MIN(v,head),MAX(v,head))) L2uh++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)INPUT_PARAM[j];
          CHANGE_STAT[j] += ((L2uh + echange == deg)
          - (L2uh == deg));
        }
      }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);

  /* *** don't forget tail -> head */    
  FOR_EACH_TOGGLE(i) {
    L2th=0;
    echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 1;
    STEP_THROUGH_OUTEDGES(head, e, u) {
      if (IS_OUTEDGE(MIN(u,tail), MAX(u,tail))){
        L2th++;
        L2tu=0;
        L2uh=0;
        STEP_THROUGH_OUTEDGES(u, f, v) {
          if(IS_OUTEDGE(MIN(v,head),MAX(v,head))) L2uh++;
          if(IS_OUTEDGE(MIN(v,tail),MAX(v,tail))) L2tu++;
        }
        STEP_THROUGH_INEDGES(u, f, v) {
          if(IS_OUTEDGE(MIN(v,head),MAX(v,head))) L2uh++;
          if(IS_OUTEDGE(MIN(v,tail),MAX(v,tail))) L2tu++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)INPUT_PARAM[j];
          CHANGE_STAT[j] -= ((L2tu + echange == deg)
          - (L2tu == deg));
          CHANGE_STAT[j] -= ((L2uh + echange == deg)
          - (L2uh == deg));
        }
      }
    }
    STEP_THROUGH_INEDGES(head, e, u) {
      if (IS_OUTEDGE(MIN(u,tail), MAX(u,tail))){
        L2th++;
        L2tu=0;
        L2uh=0;
        STEP_THROUGH_OUTEDGES(u, f, v) {
          if(IS_OUTEDGE(MIN(v,head),MAX(v,head))) L2uh++;
          if(IS_OUTEDGE(MIN(v,tail),MAX(v,tail))) L2tu++;
        }
        STEP_THROUGH_INEDGES(u, f, v) {
          if(IS_OUTEDGE(MIN(v,head),MAX(v,head))) L2uh++;
          if(IS_OUTEDGE(MIN(v,tail),MAX(v,tail))) L2tu++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)INPUT_PARAM[j];
          CHANGE_STAT[j] -= ((L2tu + echange == deg)
          - (L2tu == deg));
          CHANGE_STAT[j] -= ((L2uh + echange == deg)
          - (L2uh == deg));
        }
      }
    }
    for(j = 0; j < N_CHANGE_STATS; j++){
      deg = (Vertex)INPUT_PARAM[j];
/*      CHANGE_STAT[j] += echange*((L2th == deg) - (0 == deg)); */
      CHANGE_STAT[j] -= echange*(L2th == deg);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }  
  UNDO_PREVIOUS_TOGGLES(i);
}

/********************  changestats:  O    ***********/

// A macro indicating whether x is in [from,to)
#define FROM_TO(x, from, to) ((x)>=(from) && (x)<(to))

/*****************
 changestat: d_odegrange
*****************/
D_CHANGESTAT_FN(d_odegrange) { 
  int i, j, echange;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    Vertex tail, head;
    echange=(EdgetreeSearch(tail=TAIL(i), head=HEAD(i), nwp->outedges)==0)? 1:-1;
    Vertex tailodeg = OUT_DEG[tail];
    for(j = 0; j < N_CHANGE_STATS; j++) {
      Vertex from = INPUT_PARAM[2*j], to = INPUT_PARAM[2*j+1];
      CHANGE_STAT[j] += FROM_TO(tailodeg + echange, from, to) - FROM_TO(tailodeg, from, to);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}
 
/*****************
 changestat: d_odegrange_by_attr
*****************/
D_CHANGESTAT_FN(d_odegrange_by_attr) { 
  /* The inputparams are assumed to be set up as follows:
  The first 3*nstats values are in triples:  (from, to, attrvalue)
  The values following the first 2*nstats values are the nodal attributes.
  */
  int i, j;
  Vertex *od;
  
  od=OUT_DEG;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    Vertex tail, head;
    int echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1:1;
    Vertex tailodeg = od[tail];
    int tailattr = INPUT_PARAM[3*N_CHANGE_STATS + tail - 1]; 
    for(j = 0; j < N_CHANGE_STATS; j++){
      Vertex from = INPUT_PARAM[3*j], to = INPUT_PARAM[3*j + 1];
      int testattr = INPUT_PARAM[3*j + 2]; 
      if (tailattr == testattr)  /* we have tail attr match */
        CHANGE_STAT[j] += FROM_TO(tailodeg + echange, from, to) - FROM_TO(tailodeg, from, to);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_odegrange_w_homophily
*****************/
D_CHANGESTAT_FN(d_odegrange_w_homophily) { 
  /*  The inputparams are assumed to be set up as follows:
  The first 2*nstats values are the values of odegrange
  The values following the first 2*nstats values are the nodal attributes.
  */
  int i, j;
  double *nodeattr;
  Edge e;

  nodeattr = mtp->inputparams + N_CHANGE_STATS*2 - 1;  

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    Vertex tail=TAIL(i), head=HEAD(i);
    int tailattr = nodeattr[tail], headattr = nodeattr[head];
    if (tailattr == headattr) { /* They match; otherwise don't bother */
      int echange = IS_OUTEDGE(tail, head) ? -1:1;
      Vertex tailodeg=0, v;
      STEP_THROUGH_OUTEDGES(tail, e, v) { tailodeg += (nodeattr[v]==tailattr); }
      for(j = 0; j < N_CHANGE_STATS; j++) {
        Vertex from = INPUT_PARAM[2*j], to = INPUT_PARAM[2*j+1];
        CHANGE_STAT[j] += FROM_TO(tailodeg + echange, from, to) - FROM_TO(tailodeg, from, to);
      }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}                                        

#undef FROM_TO

/*****************
 changestat: d_odegree
*****************/
D_CHANGESTAT_FN(d_odegree) { 
  int i, j;
  
  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    Vertex tail, head;
    int echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 1;
    Vertex taild = OUT_DEG[tail];
    
    for(j=0; j < N_CHANGE_STATS; j++) {
      Vertex deg = INPUT_PARAM[j];
      CHANGE_STAT[j] += (taild + echange == deg) - (taild == deg);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_odegree_by_attr
*****************/
D_CHANGESTAT_FN(d_odegree_by_attr) { 
  /* The inputparams are assumed to be set up as follows:
  The first 2*nstats values are in pairs:  (degree, attrvalue)
  The values following the first 2*nstats values are the nodal attributes.
  */
  int i, j, echange, tailattr, testattr;
  Vertex tail, taildeg, d, *od;
  
  od=OUT_DEG;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    echange=(EdgetreeSearch(tail=TAIL(i), HEAD(i), nwp->outedges)==0)? 1:-1;
    taildeg = od[tail];
    tailattr = INPUT_PARAM[2*N_CHANGE_STATS + tail - 1]; 
    for(j = 0; j < N_CHANGE_STATS; j++) {
      d = (Vertex)INPUT_PARAM[2*j];
      testattr = INPUT_PARAM[2*j + 1]; 
      if (tailattr == testattr) { /* we have tail attr match */
        CHANGE_STAT[j] += (taildeg + echange == d) - (taildeg == d);
      }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_odegree_w_homophily
*****************/
D_CHANGESTAT_FN(d_odegree_w_homophily) { 
  /*  The inputparams are assumed to be set up as follows:
  The first nstats values are the values of degree
  The values following the first nstats values are the nodal attributes.
  */
  int i, j;
  double *nodeattr;
  Edge e;

  nodeattr = mtp->inputparams + N_CHANGE_STATS - 1;  

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    Vertex tail=TAIL(i), head=HEAD(i);
    int tailattr = nodeattr[tail], headattr = nodeattr[head];
    if (tailattr == headattr) { /* They match; otherwise don't bother */
      int echange=(EdgetreeSearch(tail, head, nwp->outedges)==0)? 1:-1;
      Vertex taildeg=0, tmp;
      STEP_THROUGH_OUTEDGES(tail, e, tmp){
        taildeg += (nodeattr[tmp]==tailattr);
      }
/*      for(e = EdgetreeMinimum(nwp->inedges, tail); */
/*      (tmp = nwp->inedges[e].value) != 0; */
/*      e = EdgetreeSuccessor(nwp->inedges, e)) { */
/*        taildeg += (nodeattr[tmp]==tailattr); */
/*      } */
      for(j = 0; j < N_CHANGE_STATS; j++) {
        Vertex deg = INPUT_PARAM[j];
        CHANGE_STAT[j] += (taildeg + echange == deg) - (taildeg == deg);
      }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_opentriad
*****************/
D_CHANGESTAT_FN(d_opentriad) { 
  int i;
  
  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    Vertex tail = TAIL(i), head = HEAD(i), node3;
    Edge change = 0, e;
    /* edgeflag is 1 if edge exists and will disappear
       edgeflag is 0 if edge DNE and will appear */
    unsigned int edgeflag = IS_OUTEDGE(tail, head);

    // -3 * triangles

    STEP_THROUGH_OUTEDGES(head, e, node3) { /* step through outedges of head */
      change += IS_UNDIRECTED_EDGE(node3,tail);
    }
    STEP_THROUGH_INEDGES(head, e, node3) { /* step through inedges of head */
      change += IS_UNDIRECTED_EDGE(node3,tail);
    }
    CHANGE_STAT[0] += change * (edgeflag ? 3.0 : -3.0);
    

    // +1 * 2-stars
    
    Vertex taild = OUT_DEG[tail] + IN_DEG[tail] - edgeflag; 
    Vertex headd = OUT_DEG[head] + IN_DEG[head] - edgeflag;
    change = taild + headd; 
    CHANGE_STAT[0] += (edgeflag ?  -change : change); 

    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_ostar
*****************/
D_CHANGESTAT_FN(d_ostar) { 
  double change, headd=0.0;
  int edgeflag, i, j, kmo;
  Edge e;
  Vertex tail, head, node3;
  int ninputs, nstats;
  double headattr;
  
  ninputs = (int)N_INPUT_PARAMS;
  nstats  = (int)N_CHANGE_STATS;
  
  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  if(ninputs>nstats){
    /* match on attributes */
    for (i=0; i < ntoggles; i++) {
      /* edgeflag is 1 if edge exists and will disappear
      edgeflag is 0 if edge DNE and will appear */
      edgeflag = IS_OUTEDGE(tail = TAIL(i), head = HEAD(i));
      headattr = INPUT_ATTRIB[head-1];
      if(headattr == INPUT_ATTRIB[tail-1]){
        headd = - edgeflag;
        for(e = EdgetreeMinimum(nwp->outedges, tail);
	      (node3 = nwp->outedges[e].value) != 0;
	      e = EdgetreeSuccessor(nwp->outedges, e)) { /* step through outedges of head */
          if(headattr == INPUT_ATTRIB[node3-1]){++headd;}
        }
        for(j=0; j < N_CHANGE_STATS; j++) {
          kmo = ((int)INPUT_PARAM[j]) - 1;
          change = CHOOSE(headd, kmo); 
          CHANGE_STAT[j] += (edgeflag ? - change : change); 
        }
      }
    TOGGLE_IF_MORE_TO_COME(i);
    }
  }else{
    for (i=0; i < ntoggles; i++) {
      /* edgeflag is 1 if edge exists and will disappear
      edgeflag is 0 if edge DNE and will appear */
      edgeflag = IS_OUTEDGE(tail = TAIL(i), head = HEAD(i));
      headd = OUT_DEG[tail] - edgeflag;      
      for(j=0; j < N_CHANGE_STATS; j++) {
        kmo = ((int)INPUT_PARAM[j]) - 1;
        change = CHOOSE(headd, kmo); 
        CHANGE_STAT[j] += (edgeflag ? - change : change); 
      }
    TOGGLE_IF_MORE_TO_COME(i);
    }
  }
  
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_odegreepopularity
*****************/
D_CHANGESTAT_FN(d_odegreepopularity) { 
  int i, edgeflag;
  double change;
  Vertex head, tail, deg=0;
  
  /* *** don't forget tail -> head */    
  change = 0.0;
  FOR_EACH_TOGGLE(i) {
    tail=TAIL(i);
    head=HEAD(i);
    edgeflag = IS_OUTEDGE(tail, head); /* either 0 or 1 */
    deg = (double)(OUT_DEG[tail]);
    if(edgeflag){
      change -= sqrt(deg);
      change += (deg-1.0)*(sqrt(deg-1.0)-sqrt(deg));
    }else{
      change += sqrt(deg+1.0);
      change += deg*(sqrt(deg+1.0)-sqrt(deg));
    }
    TOGGLE_IF_MORE_TO_COME(i); 
  }
  CHANGE_STAT[0]=change; 
  UNDO_PREVIOUS_TOGGLES(i);
}

/********************  changestats:  R    ***********/
/*****************
 changestat: d_receiver
*****************/
D_CHANGESTAT_FN(d_receiver) { 
  int i, j, echange;
  Vertex tail, head, deg;
  
  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {      
    echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 1;
    j=0;
    deg = (Vertex)INPUT_PARAM[j];
    while((deg != head) && (j < (N_CHANGE_STATS-1))){
      j++;
      deg = (Vertex)INPUT_PARAM[j];
    }
    if(deg==head){CHANGE_STAT[j] += echange;}
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/********************  changestats:  S    ***********/
/*****************
 changestat: d_sender
*****************/
D_CHANGESTAT_FN(d_sender) { 
  int i, j, echange;
  Vertex tail, head, deg;
  
  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 1;
    j=0;
    deg = (Vertex)INPUT_PARAM[j];
    while((deg != tail) && (j < (N_CHANGE_STATS-1))){
      j++;
      deg = (Vertex)INPUT_PARAM[j];
    }
    if(deg==tail){CHANGE_STAT[j] += echange;}
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_simmelian
*****************/
D_CHANGESTAT_FN(d_simmelian) { 
  Edge e;
  Vertex tail, head, change, node3;
  int edgeflag, i;
  
  /* *** don't forget tail -> head */    
 CHANGE_STAT[0] = 0.0;
 FOR_EACH_TOGGLE(i) 
 {
  edgeflag = IS_OUTEDGE(tail = TAIL(i), head = HEAD(i));
   
  if(EdgetreeSearch(head, tail, nwp->outedges) != 0){
   change = 0;
   
   for(e = EdgetreeMinimum(nwp->outedges, head);
       (node3 = nwp->outedges[e].value) != 0;
       e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of head */
   {
     if (node3 != tail
      && EdgetreeSearch(node3, tail, nwp->outedges) != 0 
      && EdgetreeSearch(tail, node3, nwp->outedges) != 0 
      && EdgetreeSearch(node3, head, nwp->outedges) != 0 
        ){++change;}
   }
      
   CHANGE_STAT[0] += edgeflag ? -(double)change : (double)change;
   }
   
   TOGGLE_IF_MORE_TO_COME(i);
  }
  
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_simmelianties
*****************/
D_CHANGESTAT_FN(d_simmelianties) { 
  Edge e, e2;
  Vertex tail, head, change, node3, node4, first, htflag;
  int edgeflag, i;
  
  /* *** don't forget tail -> head */    
 CHANGE_STAT[0] = 0.0;
 FOR_EACH_TOGGLE(i) {
   edgeflag = IS_OUTEDGE(tail = TAIL(i), head = HEAD(i));

   if(IS_OUTEDGE(head, tail)){
     change = htflag = 0;
     STEP_THROUGH_OUTEDGES(head, e, node3) { /* step through outedges of head */
       if (node3 != tail
       && IS_OUTEDGE(node3, tail) && IS_OUTEDGE(tail, node3) && IS_OUTEDGE(node3, head)){
         htflag=1; /* tail, head is itself in a simmelian triple (along with head, tail)*/
         first = 1;
         /* Find out whether (tail, node3) is in any other simmelian triple */
         STEP_THROUGH_OUTEDGES (tail, e2, node4) { /* step through outedges of tail */
           if (node4 != head  && node4 != node3 && IS_OUTEDGE(node4, tail) 
           && IS_OUTEDGE(node4, node3) && IS_OUTEDGE(node3, node4)){
             first = 0;
           }
         }
         if(first){++change;}
         first = 1;
         /* Find out whether (head, node3) is in any other simmelian triple */
         STEP_THROUGH_OUTEDGES (head, e2, node4) { /* step through outedges of head */
           if (node4 != tail  && node4 != node3 && IS_OUTEDGE(node4, head) 
           && IS_OUTEDGE(node4, node3) && IS_OUTEDGE(node3, node4)) {
             first = 0;
           }
         }
         if(first){++change;}
       }
     }
     change += htflag;
     change = 2*change; /* All changes must happen in pairs here; no tie can
                           be counted without its opposite */
     CHANGE_STAT[0] += edgeflag ? -(double)change : (double)change;
   }
   TOGGLE_IF_MORE_TO_COME(i);
 }
 UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_smalldiff
*****************/
D_CHANGESTAT_FN(d_smalldiff) { 
  Vertex tail, head;
  int i;
  
  /* *** don't forget tail -> head */    
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i) {
    tail=TAIL(i);
    head=HEAD(i);
    CHANGE_STAT[0] += (fabs(INPUT_ATTRIB[tail-1] - INPUT_ATTRIB[head-1])
    > INPUT_PARAM[0]) ? 0.0 :
    ((EdgetreeSearch(tail, head, nwp->outedges) != 0) ? -1.0 : 1.0); 
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_sociality
*****************/
D_CHANGESTAT_FN(d_sociality) { 
  int i, j, echange;
  Vertex tail, head, deg;
  int ninputs, nstats;
  double tailattr;
  
  ninputs = (int)N_INPUT_PARAMS;
  nstats  = (int)N_CHANGE_STATS;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  if(ninputs>nstats+1){
    /* match on attributes */
    FOR_EACH_TOGGLE(i) {      
      echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 1;
      tailattr = INPUT_ATTRIB[tail-1+nstats+1]; // +1 for the "guard" value between vertex IDs and attribute vector
      if(tailattr == INPUT_ATTRIB[head-1+nstats+1]){
	j=0;
	deg = (Vertex)INPUT_PARAM[j];
	while(deg != tail && j < nstats){
	  j++;
	  deg = (Vertex)INPUT_PARAM[j];
	}
	if(j < nstats){CHANGE_STAT[j] += echange;}
	j=0;
	deg = (Vertex)INPUT_PARAM[j];
	while(deg != head && j < nstats){
	  j++;
	  deg = (Vertex)INPUT_PARAM[j];
	}
	if(j < nstats){CHANGE_STAT[j] += echange;}
      }
      
      TOGGLE_IF_MORE_TO_COME(i);
    }
  }else{
    /* *** don't forget tail -> head */    
    FOR_EACH_TOGGLE(i) {      
      echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 1;
      j=0;
      deg = (Vertex)INPUT_PARAM[j];
      while(deg != tail && j < nstats){
	j++;
	deg = (Vertex)INPUT_PARAM[j];
      }
      if(j < nstats){CHANGE_STAT[j] += echange;}
      j=0;
      deg = (Vertex)INPUT_PARAM[j];
      while(deg != head && j < nstats){
	j++;
	deg = (Vertex)INPUT_PARAM[j];
      }
      if(j < nstats){CHANGE_STAT[j] += echange;}
      
      TOGGLE_IF_MORE_TO_COME(i);
    }
  }
  
  UNDO_PREVIOUS_TOGGLES(i);
}

/********************  changestats:  T    ***********/
/*****************
 changestat: d_tdsp
*****************/
D_CHANGESTAT_FN(d_tdsp) {
  Edge e, f;
  int i, j, echange, L2tu, L2uh;
  Vertex deg, tail, head, u, v;
  
  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i){
    tail = TAIL(i); head=HEAD(i);
    echange = 1-2*IS_OUTEDGE(tail,head);
    /* step through outedges of head */
    for(e = MIN_OUTEDGE(head); (u=OUTVAL(e))!=0; e=NEXT_OUTEDGE(e)) { 
      if (u != tail){
        L2tu=0; /* This will be # of shared partners of (tail,u) */
        /* step through inedges of u, incl. (head,u) itself */
        for(f = MIN_INEDGE(u); (v=INVAL(f))!=0; f=NEXT_INEDGE(f)) {
          if(IS_OUTEDGE(tail,v)) L2tu++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)INPUT_PARAM[j];
          CHANGE_STAT[j] += ((L2tu + echange == deg) - (L2tu == deg));
        }
      }
    }
    /* step through inedges of tail */
    for(e = MIN_INEDGE(tail); (u=INVAL(e))!=0; e=NEXT_INEDGE(e)) {
      if (u != head){
        L2uh=0; /* This will be # of shared partners of (u,head) */
        /* step through outedges of u , incl. (u,tail) itself */
        for(f = MIN_OUTEDGE(u);(v=OUTVAL(f))!=0; f=NEXT_OUTEDGE(f)){
          if(IS_OUTEDGE(v,head)) L2uh++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)INPUT_PARAM[j];
          CHANGE_STAT[j] += ((L2uh + echange == deg) - (L2uh == deg));
        }
      }
    }
    
    if (i+1 < ntoggles) TOGGLE(tail,head);  /* Toggle this edge if more to come */
  }
  
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_tesp
*****************/
D_CHANGESTAT_FN(d_tesp) { 
  Edge e, f;
  int i, j, echange;
  int L2th, L2tu, L2uh;
  Vertex deg;
  Vertex tail, head, u, v;
  
  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i){      
    L2th=0;
    echange = (EdgetreeSearch(tail=TAIL(i), head=HEAD(i), nwp->outedges) == 0) ? 1 : -1;
    /* step through outedges of head */
    for(e = EdgetreeMinimum(nwp->outedges, head); 
    (u = nwp->outedges[e].value) != 0; e = EdgetreeSuccessor(nwp->outedges, e)) {
      if (EdgetreeSearch(tail, u, nwp->outedges) != 0){
        L2tu=0;
        /* step through inedges of u */
        for(f = EdgetreeMinimum(nwp->inedges, u); (v = nwp->inedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->inedges, f)){
          if(EdgetreeSearch(tail,v,nwp->outedges)!= 0) L2tu++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)INPUT_PARAM[j];
          CHANGE_STAT[j] += ((L2tu + echange == deg) - (L2tu == deg));
        }
      }
    }
    /* step through inedges of head */
    for (e = EdgetreeMinimum(nwp->inedges, head); (u = nwp->inedges[e].value) != 0;
    e = EdgetreeSuccessor(nwp->inedges, e)){
      if (EdgetreeSearch(tail, u, nwp->outedges) != 0){
        L2th++;
      }
      if (EdgetreeSearch(u, tail, nwp->outedges) != 0){
        L2uh=0;
        /* step through outedges of u */
        for(f = EdgetreeMinimum(nwp->outedges, u); (v = nwp->outedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->outedges, f)){
          if(EdgetreeSearch(v, head, nwp->outedges)!= 0) L2uh++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)INPUT_PARAM[j];
          CHANGE_STAT[j] += ((L2uh + echange == deg) - (L2uh == deg));
        }
      }
    }
    for(j = 0; j < N_CHANGE_STATS; j++){
      deg = (Vertex)INPUT_PARAM[j];
      CHANGE_STAT[j] += echange*(L2th == deg);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }  
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_threepath
*****************/
D_CHANGESTAT_FN(d_threepath) { 
  int i, j, k, edgeflag, change, dchange[4];
  Edge e;
  Vertex tail, head, node3;
  /* The four values of dchange represent the four different types of
     directed threepaths oriented so that the middle step is always 
     "right" (R).  In order:   RRR, RRL, LRR, LRL 
                         i.e., >>>  >><  <>>  <><  */


  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    edgeflag = IS_OUTEDGE(tail = TAIL(i), head = HEAD(i));
    /* Step A: Count threepaths in which tail->head is the middle edge */
    dchange[0] = IN_DEG[tail] * OUT_DEG[head]; /* R then R; may count head->tail->head->tail */
    dchange[1] = IN_DEG[tail] * (IN_DEG[head]-edgeflag); /* R then L */
    dchange[2] = (OUT_DEG[tail]-edgeflag) * OUT_DEG[head]; /* L then R */
    dchange[3] = (OUT_DEG[tail]-edgeflag) * (IN_DEG[head]-edgeflag); /* L then L */
    /* Step B: Count threepaths where tail is one endpoint */
    STEP_THROUGH_OUTEDGES(head, e, node3) { /* tail->head->node3-x  which means -RL */
      dchange[1] += IN_DEG[node3]-1;    /* RRL; subtract 1 for head itself */
      dchange[0] += OUT_DEG[node3]; /* RRR; possibly counted tail->head->tail->head */
    }
    STEP_THROUGH_INEDGES(head, e, node3) { /* x-node3->head<-tail  which means -RL*/
      if (node3 != tail) {
        dchange[3] += OUT_DEG[node3]-1;  /* LRL; subtract 1 for head itself */
        dchange[1] += IN_DEG[node3];     /* RRL */
      }
    }
    /* Step C: Count threepaths where head is one endpoint */
    STEP_THROUGH_INEDGES(tail, e, node3) { /* x-node3->tail->head  which means -RR */
      dchange[2] += OUT_DEG[node3]-1;  /* LRR; subtract 1 for tail itself */
      dchange[0] += IN_DEG[node3]; /* RRR; possibly counted tail->head->tail->head */
    }
    STEP_THROUGH_OUTEDGES(tail, e, node3) { /* head<-tail->node3-x  which means LR- */
      if (node3 != head) {
        dchange[3] += IN_DEG[node3]-1; /* LRL; subtract 1 for tail itself */
        dchange[2] += OUT_DEG[node3];  /* LRR */
      }
    }
    /* Finally, correct for overcounted head->tail->head->tail and tail->head->tail->head */
    if (DIRECTED) {
      dchange[0] -= IS_INEDGE(tail, head) * (1 + 2 * edgeflag);
      /* head->tail->head->tail is counted in A whenever IS_INEDGE(tail,head) but 
         TT->head->tail->head is only counted in B and C when also IS_OUTEDGE(tail, head) */
      for (j = 0; j < N_INPUT_PARAMS; j++) {
        k = (int) INPUT_PARAM[j];
        CHANGE_STAT[j] += (edgeflag ? -dchange[k-1] : dchange[k-1]);
      }
    }
    else { /* Undirected case; don't need head->tail->head->tail correction */
      change = dchange[0] + dchange[1] + dchange[2] + dchange[3];
      CHANGE_STAT[0] += (edgeflag ? -change : change); 
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_tnsp
*****************/
D_CHANGESTAT_FN(d_tnsp) { 
  Edge e, f;
  int i, j, echange;
  int L2th, L2tu, L2uh;
  Vertex deg;
  Vertex tail, head, u, v;
  
  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i){
    tail = TAIL(i); head=HEAD(i);
    echange = 1-2*IS_OUTEDGE(tail,head);
    /* step through outedges of head */
    for(e = MIN_OUTEDGE(head); (u=OUTVAL(e))!=0; e=NEXT_OUTEDGE(e)) { 
      if (u != tail){
        L2tu=0; /* This will be # of shared partners of (tail,u) */
        /* step through inedges of u, incl. (head,u) itself */
        for(f = MIN_INEDGE(u); (v=INVAL(f))!=0; f=NEXT_INEDGE(f)) {
          if(IS_OUTEDGE(tail,v)) L2tu++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)INPUT_PARAM[j];
          CHANGE_STAT[j] += ((L2tu + echange == deg) - (L2tu == deg));
        }
      }
    }
    /* step through inedges of tail */
    for(e = MIN_INEDGE(tail); (u=INVAL(e))!=0; e=NEXT_INEDGE(e)) {
      if (u != head){
        L2uh=0; /* This will be # of shared partners of (u,head) */
        /* step through outedges of u , incl. (u,tail) itself */
        for(f = MIN_OUTEDGE(u);(v=OUTVAL(f))!=0; f=NEXT_OUTEDGE(f)){
          if(IS_OUTEDGE(v,head)) L2uh++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)INPUT_PARAM[j];
          CHANGE_STAT[j] += ((L2uh + echange == deg) - (L2uh == deg));
        }
      }
    }
    
    if (i+1 < ntoggles) TOGGLE(tail, head);  /* Toggle this edge if more to come */
  }
  
  UNDO_PREVIOUS_TOGGLES(i);

    /* *** don't forget tail -> head */    
    FOR_EACH_TOGGLE(i){      
    L2th=0;
    echange = (EdgetreeSearch(tail=TAIL(i), head=HEAD(i), nwp->outedges) == 0) ? 1 : -1;
    /* step through outedges of head */
    for(e = EdgetreeMinimum(nwp->outedges, head); 
    (u = nwp->outedges[e].value) != 0; e = EdgetreeSuccessor(nwp->outedges, e)) {
      if (EdgetreeSearch(tail, u, nwp->outedges) != 0){
        L2tu=0;
        /* step through inedges of u */
        for(f = EdgetreeMinimum(nwp->inedges, u); (v = nwp->inedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->inedges, f)){
          if(EdgetreeSearch(tail,v,nwp->outedges)!= 0) L2tu++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)INPUT_PARAM[j];
          CHANGE_STAT[j] -= ((L2tu + echange == deg) - (L2tu == deg));
        }
      }
    }
    /* step through inedges of head */
    for (e = EdgetreeMinimum(nwp->inedges, head); (u = nwp->inedges[e].value) != 0;
    e = EdgetreeSuccessor(nwp->inedges, e)){
      if (EdgetreeSearch(tail, u, nwp->outedges) != 0){
        L2th++;
      }
      if (EdgetreeSearch(u, tail, nwp->outedges) != 0){
        L2uh=0;
        /* step through outedges of u */
        for(f = EdgetreeMinimum(nwp->outedges, u); (v = nwp->outedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->outedges, f)){
          if(EdgetreeSearch(v, head, nwp->outedges)!= 0) L2uh++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)INPUT_PARAM[j];
          CHANGE_STAT[j] -= ((L2uh + echange == deg) - (L2uh == deg));
        }
      }
    }
    for(j = 0; j < N_CHANGE_STATS; j++){
      deg = (Vertex)INPUT_PARAM[j];
      CHANGE_STAT[j] -= echange*(L2th == deg);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }  
  UNDO_PREVIOUS_TOGGLES(i);

}

/*****************
 changestat: d_transitive
*****************/
D_CHANGESTAT_FN(d_transitive) { 
  Edge e;
  Vertex tail, head, node2;
  double change;
  int edgeflag, i;
  
  /* *** don't forget tail -> head */    
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i) {
    edgeflag = IS_OUTEDGE(tail = TAIL(i), head = HEAD(i));
    change = 0.0; /* change should become the number of transitive triples
                     a->b, b->c, a->c in which tail->head is found  */
    
    STEP_THROUGH_OUTEDGES(head, e, node2) { /* step through outedges of head */
      if (tail != node2 && IS_OUTEDGE(tail, node2)){
        change = change + 1.0; /* Here we have tail->head, head->node2, tail->node2 */
      }
    }
    STEP_THROUGH_INEDGES(head, e, node2) { /* step through inedges of head */
      if (tail != node2) {
        change = change + IS_OUTEDGE(tail, node2) + IS_OUTEDGE(node2, tail);
        /* Here we have tail->head and node2->head, with either node2->tail or tail->node2 */
      }
    }
//    STEP_THROUGH_INEDGES(tail, e, node2) { /* step through inedges of tail */
//      if (node2 != head){
//        if (!IS_OUTEDGE(node2, head)){
//          change = change - 1.0;
//        }
//      }
//    }
    CHANGE_STAT[0] += edgeflag ? -change : change;
//  Rprintf("tail %d head %d edgeflag %d change %f C_S[0]=%f\n", tail, head, edgeflag, change,CHANGE_STAT[0]); 
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

D_CHANGESTAT_FN(d_transitiveties) { 
  Edge e, f;
  int i, echange, ochange;
  int L2th, L2tu, L2uh;
  Vertex tail, head, u, v;
  double cumchange;
  double tailattr;
  
  CHANGE_STAT[0] = 0.0;
  
  /* *** don't forget tail -> head */    
  FOR_EACH_TOGGLE(i){      
    cumchange=0.0;
    L2th=0;
    ochange = (EdgetreeSearch(tail=TAIL(i), head=HEAD(i), nwp->outedges) == 0) ? 0 : -1;
    echange = 2*ochange + 1;
    if(N_INPUT_PARAMS>0){ /* match on attributes */
      tailattr = INPUT_ATTRIB[tail-1];
      if(tailattr == INPUT_ATTRIB[head-1]){
       /* step through outedges of head  */
       for(e = EdgetreeMinimum(nwp->outedges, head);
	   (u = nwp->outedges[e].value) != 0;
	   e = EdgetreeSuccessor(nwp->outedges, e)){
         if (EdgetreeSearch(tail, u, nwp->outedges) != 0 && (tailattr == INPUT_ATTRIB[u-1])){
	   L2tu=ochange;
	   /* step through inedges of u */
	   for(f = EdgetreeMinimum(nwp->inedges, u); 
	       (v = nwp->inedges[f].value) != 0;
	       f = EdgetreeSuccessor(nwp->inedges, f)){
	     if(EdgetreeSearch(tail,v,nwp->outedges)!= 0 && (tailattr == INPUT_ATTRIB[v-1])){
	       L2tu++;
	       if(L2tu>0) {break;}
	     }
	   }
	   cumchange += (L2tu==0);
         }
       }
       /* step through inedges of head */
       
       for(e = EdgetreeMinimum(nwp->inedges, head);
	   (u = nwp->inedges[e].value) != 0;
	   e = EdgetreeSuccessor(nwp->inedges, e)){
         if (EdgetreeSearch(tail, u, nwp->outedges) != 0 && (tailattr == INPUT_ATTRIB[u-1])){
	   L2th++;
         }
         if (EdgetreeSearch(u, tail, nwp->outedges) != 0 && (tailattr == INPUT_ATTRIB[u-1])){
	   L2uh=ochange;
	   /* step through outedges of u */
	   for(f = EdgetreeMinimum(nwp->outedges, u);
	       (v = nwp->outedges[f].value) != 0;
	       f = EdgetreeSuccessor(nwp->outedges, f)){
	     if(EdgetreeSearch(v,head,nwp->outedges)!= 0 && (tailattr == INPUT_ATTRIB[v-1])){
	       L2uh++;
	       if(L2uh>0) {break;}
	     }
	   }
	   cumchange += (L2uh==0) ;
         }
       }}
      }else{ /* no attributes */
    /* step through outedges of head  */
    for(e = EdgetreeMinimum(nwp->outedges, head);
	(u = nwp->outedges[e].value) != 0;
	e = EdgetreeSuccessor(nwp->outedges, e)){
      if (EdgetreeSearch(tail, u, nwp->outedges) != 0){
	L2tu=ochange;
	/* step through inedges of u */
	for(f = EdgetreeMinimum(nwp->inedges, u); 
	    (v = nwp->inedges[f].value) != 0;
	    f = EdgetreeSuccessor(nwp->inedges, f)){
	  if(EdgetreeSearch(tail,v,nwp->outedges)!= 0){
	    L2tu++;
	    if(L2tu>0) {break;}
	  }
	}
	cumchange += (L2tu==0);
      }
    }
    /* step through inedges of head */
    
    for(e = EdgetreeMinimum(nwp->inedges, head);
	(u = nwp->inedges[e].value) != 0;
	e = EdgetreeSuccessor(nwp->inedges, e)){
      if (EdgetreeSearch(tail, u, nwp->outedges) != 0){
	L2th++;
      }
      if (EdgetreeSearch(u, tail, nwp->outedges) != 0){
	L2uh=ochange;
	/* step through outedges of u */
	for(f = EdgetreeMinimum(nwp->outedges, u);
	    (v = nwp->outedges[f].value) != 0;
	    f = EdgetreeSuccessor(nwp->outedges, f)){
	  if(EdgetreeSearch(v,head,nwp->outedges)!= 0){
	    L2uh++;
	    if(L2uh>0) {break;}
	  }
	}
	cumchange += (L2uh==0) ;
      }
    }
    }
    
    cumchange += (L2th>0) ;
//  Rprintf("L2th %d echange %d cumchange %f tail %d head %d\n", L2th, echange, cumchange,tail,head);
    cumchange  = echange*cumchange;
    (CHANGE_STAT[0]) += cumchange;
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}
///*****************
// globalstat: s_transitiveties
//*****************/
//S_CHANGESTAT_FN(s_transitiveties) { 
//  Edge e1, e2;
//  Vertex tail, head, change, node3;
//  double tailattr;
//  int hnottrans;
//  
//  /* *** don't forget tail -> head */    
//  change=0;
//  if(N_INPUT_PARAMS > 0){ /* match on attributes */
//    for (tail=1; tail <= N_NODES; tail++) { 
//      tailattr = INPUT_ATTRIB[tail-1];
//      STEP_THROUGH_OUTEDGES(tail, e1, head) {
//        if(tailattr == INPUT_ATTRIB[head-1]) {
//          hnottrans=1;
//          STEP_THROUGH_INEDGES(head, e2, node3) { 
//            if(hnottrans && IS_INEDGE(node3, tail) && (tailattr == INPUT_ATTRIB[node3-1])){ /* tail -> head base forms transitive */
//              hnottrans=0;
//              change++;
//            }
//          }
//        }
//      }
//    }
//  }else{
//    /* *** don't forget tail -> head */    
//    for (tail=1; tail <= N_NODES; tail++) { 
//      STEP_THROUGH_OUTEDGES(tail, e1, head) {
//        hnottrans=1;
//        STEP_THROUGH_INEDGES(head, e2, node3) { 
//          if(hnottrans && IS_INEDGE(node3, tail)){ /* tail -> head base forms transitive */
//            hnottrans=0;
//            change++;
//          }
//        }
//      }
//    }
//  }
//  CHANGE_STAT[0] = change;
//}

D_CHANGESTAT_FN(d_cyclicalties) { 
  Edge e, f;
  int i, echange, ochange;
  int L2th, L2tu, L2uh;
  Vertex tail, head, u, v;
  double cumchange;
  double tailattr;
  
  CHANGE_STAT[0] = 0.0;
  
  /* *** don't forget tail -> head */    
  FOR_EACH_TOGGLE(i){      
    cumchange=0.0;
    L2th=0;
    ochange = (EdgetreeSearch(tail=TAIL(i), head=HEAD(i), nwp->outedges) == 0) ? 0 : -1;
    echange = 2*ochange + 1;
    if(N_INPUT_PARAMS>0){ /* match on attributes */
      tailattr = INPUT_ATTRIB[tail-1];
      if(tailattr == INPUT_ATTRIB[head-1]){
       /* step through outedges of head  */
       for(e = EdgetreeMinimum(nwp->outedges, head);
	   (u = nwp->outedges[e].value) != 0;
	   e = EdgetreeSuccessor(nwp->outedges, e)){
         if (EdgetreeSearch(tail, u, nwp->inedges) != 0 && (tailattr == INPUT_ATTRIB[u-1])){
	   L2tu=ochange;
	   /* step through inedges of u */
	   for(f = EdgetreeMinimum(nwp->inedges, u); 
	       (v = nwp->inedges[f].value) != 0;
	       f = EdgetreeSuccessor(nwp->inedges, f)){
	     if(EdgetreeSearch(tail,v,nwp->outedges)!= 0 && (tailattr == INPUT_ATTRIB[v-1])){
	       L2tu++;
	       if(L2tu>0) {break;}
	     }
	   }
	   cumchange += (L2tu==0);
         }
       }
       /* step through inedges of head */
       
       for(e = EdgetreeMinimum(nwp->outedges, head);
	   (u = nwp->outedges[e].value) != 0;
	   e = EdgetreeSuccessor(nwp->outedges, e)){
         if (EdgetreeSearch(tail, u, nwp->inedges) != 0 && (tailattr == INPUT_ATTRIB[u-1])){
	   L2th++;
         }
         if (EdgetreeSearch(u, tail, nwp->outedges) != 0 && (tailattr == INPUT_ATTRIB[u-1])){
	   L2uh=ochange;
	   /* step through outedges of u */
	   for(f = EdgetreeMinimum(nwp->outedges, u);
	       (v = nwp->outedges[f].value) != 0;
	       f = EdgetreeSuccessor(nwp->outedges, f)){
	     if(EdgetreeSearch(v,head,nwp->outedges)!= 0 && (tailattr == INPUT_ATTRIB[v-1])){
	       L2uh++;
	       if(L2uh>0) {break;}
	     }
	   }
	   cumchange += (L2uh==0) ;
         }
       }}
      }else{ /* no attributes */
    /* step through outedges of head  */
    for(e = EdgetreeMinimum(nwp->outedges, head);
	(u = nwp->outedges[e].value) != 0;
	e = EdgetreeSuccessor(nwp->outedges, e)){
      if (EdgetreeSearch(tail, u, nwp->inedges) != 0){
	L2tu=ochange;
	/* step through inedges of u */
	for(f = EdgetreeMinimum(nwp->inedges, u); 
	    (v = nwp->inedges[f].value) != 0;
	    f = EdgetreeSuccessor(nwp->inedges, f)){
	  if(EdgetreeSearch(tail,v,nwp->outedges)!= 0){
	    L2tu++;
	    if(L2tu>0) {break;}
	  }
	}
	cumchange += (L2tu==0);
      }
    }
    /* step through outedges of head */
    
    for(e = EdgetreeMinimum(nwp->outedges, head);
	(u = nwp->outedges[e].value) != 0;
	e = EdgetreeSuccessor(nwp->outedges, e)){
      if (EdgetreeSearch(tail, u, nwp->inedges) != 0){
	L2th++;
      }
      if (EdgetreeSearch(u, tail, nwp->outedges) != 0){
	L2uh=ochange;
	/* step through outedges of u */
	for(f = EdgetreeMinimum(nwp->outedges, u);
	    (v = nwp->outedges[f].value) != 0;
	    f = EdgetreeSuccessor(nwp->outedges, f)){
	  if(EdgetreeSearch(v,head,nwp->outedges)!= 0){
	    L2uh++;
	    if(L2uh>0) {break;}
	  }
	}
	cumchange += (L2uh==0) ;
      }
    }
    }
    
    cumchange += (L2th>0) ;
//  Rprintf("L2th %d echange %d cumchange %f tail %d head %d\n", L2th, echange, cumchange,tail,head);
    cumchange  = echange*cumchange;
    (CHANGE_STAT[0]) += cumchange;
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_triadcensus
*****************/
D_CHANGESTAT_FN(d_triadcensus) { 
  int i, j, edgeflag, a, b, c, d, e, edgecount, t300, 
  t210, t120C, t120U, t120D, t201, t030C, t030T, t111U, 
  t111D, t021C, t021U, t021D, t102, t012, t003;
  Vertex triadtype, node3, tail, head;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  if (DIRECTED) {
    /* directed version */
    FOR_EACH_TOGGLE(i) {      
      edgeflag = IS_OUTEDGE(tail = TAIL(i), head = HEAD(i));
      t300 = 0;
      t210 = 0;
      t120C = 0;  t120U = 0;   t120D = 0;  t201 = 0;
      t030C = 0;  t030T = 0;   t111U = 0;  t111D = 0;
      t021C = 0;  t021U = 0;   t021D = 0;  t102 = 0;
      t012 = 0;
      
      if ((EdgetreeMinimum(nwp->outedges, head) != 0) || 
        (EdgetreeMinimum(nwp->inedges, head) != 0) || 
        (EdgetreeMinimum(nwp->outedges, tail) != 0) ||
        (EdgetreeMinimum(nwp->inedges, tail) != 0)) {      

          /* ****** loop through node3 ****** */
          for (node3=1; node3 <= N_NODES; node3++) { 
            if (node3 != tail && node3 != head) {
              a = (EdgetreeSearch(head, tail, nwp->outedges) != 0); 
              b = (EdgetreeSearch(head, node3, nwp->outedges) != 0);
              c = (EdgetreeSearch(node3, head, nwp->outedges) != 0);
              d = (EdgetreeSearch(node3, tail, nwp->outedges) != 0);
              e = (EdgetreeSearch(tail, node3, nwp->outedges) != 0);
              edgecount = (a + b + c + d + e);
              
              switch(edgecount) {
                case 0:   /* 012 */
                ++t012;

                case 1: {  /* 021C, 021U, 021D, 102 */
                  if ((b == 1) || (d == 1))
                    ++t021C;
                  if (c == 1)
                    ++t021U;
                  if (e == 1)
                    ++t021D;
                  if (a == 1)
                    ++t102;
                }
                break;
                
                case 2: {  /* 030C, 030T, 111U, 111D */
                  if ((b + d) == 2)       
                    ++t030C;
                  if (((b + e) == 2) || ((c + d) == 2) || ((c + e) == 2))
                    ++t030T;
                  if (((a + b) == 2) || ((a + e) == 2) || ((d + e) == 2))
                    ++t111U;
                  if (((a + c) == 2) || ((a + d) == 2) || ((b + c) == 2))
                    ++t111D;
                }
                break;
            
                case 3: {   /* 120C, 120U, 120D, 201 */
                  if (a == 1) {
                    if (((b + d) == 2) || ((c + e) == 2))
                      ++t120C;
                    if ((b + e) == 2)
                      ++t120U;
                    if ((c + d) == 2)
                      ++t120D;
                    if (((b + c) == 2) || ((d + e) == 2))
                      ++t201;
                  } else {
                    if (b == 1) {
                      if (((c + d) == 2) || ((d + e) == 2))
                        ++t120C;
                      if ((c + e) == 2)
                        ++t120D;
                    } else  {
                      ++t120U;
                    }
                  } 
                }
                break;
             
                case 4:   /* 210 */
                ++t210;
                break;
             
                case 5:   /* 300 */            
                ++t300;
                break;
              }

              switch(edgecount) {  
                case 1:   /* 102, 021D, 021U, 021C */
                --t012;
                break;

                case 2: {  /* 030C, 030T, 111U, 111D */ 
                  if (((a + c) == 2) || ((a + e) == 2) || ((b + d) == 2) || 
                    ((c + e) == 2)) 
                  --t021C;
                  if (((a + d) == 2) || ((b + e) == 2))
                    --t021U;
                  if (((a + b) == 2) || ((c + d) == 2))
                    --t021D;
                  if (((b + c) == 2) || ((d + e) == 2))
                    --t102;
                } 
                break;

                case 3: {  /* 201, 120D, 120U, 120C */
                  if (a == 1) {
                    if ((c + e) == 2)       
                      --t030C;
                    if (((c + d) == 2) || ((b + e) == 2) || ((b + d) == 2))
                      --t030T;
                    if ((b + c) == 2)
                      --t111U;
                    if ((d + e) == 2)
                      --t111D;
                  } else {
                    if (b == 1) {
                      if ((c + d) == 2)
                        --t111U;
                      if (((c + e) == 2) || ((d + e) == 2))
                        --t111D;
                    } else {
                    --t111U;
                    }
                  }
                }
                break;

                case 4: {  /* 210 */
                  if (a == 1)
                  {
                    if (((b + c + e) == 3) || ((c + d + e) == 3))
                      --t120C;
                    if ((b + c + d) == 3)
                      --t120U;
                    if ((b + d + e) == 3)
                      --t120D;
                  } else {
                    if ((b + c + d + e) == 4)
                      --t201;
                  } 
                }
                break;

                case 5:   /* 300 */            
                --t210;
                break;
              }
            }
          }    /* ******  move to next node3 ******** */
        }
        else 
          t012 = t012 + (N_NODES - 2);  

        for(j = 0; j < N_CHANGE_STATS; j++) { 
          triadtype = (Vertex)INPUT_PARAM[j]; 
          
          switch(triadtype) { /* SEARCH_ON_THIS_TO_TRACK_DOWN_TRIADCENSUS_CHANGE
                                 to undo triadcensus change, change - to plus in 
                                  next two lines: */
            case 1:  t003 = -(t300+t210+t120C+t120U+t120D+t201+t030C+t030T);
            t003 = t003-(t111U+t111D+t021C+t021U+t021D+t102+t012);
            CHANGE_STAT[j] += edgeflag ? -(double)t003 : (double)t003;
            break;
            case 2:   CHANGE_STAT[j] += edgeflag ? -(double)t012 : (double)t012;
            break;
            case 3:   CHANGE_STAT[j] += edgeflag ? -(double)t102 : (double)t102;
            break;
            case 4:   CHANGE_STAT[j] += edgeflag ? -(double)t021D : (double)t021D;
            break;
            case 5:   CHANGE_STAT[j] += edgeflag ? -(double)t021U : (double)t021U;
            break;
            case 6:   CHANGE_STAT[j] += edgeflag ? -(double)t021C : (double)t021C;
            break;
            case 7:	  CHANGE_STAT[j] += edgeflag ? -(double)t111D : (double)t111D;
            break;
            case 8:	  CHANGE_STAT[j] += edgeflag ? -(double)t111U : (double)t111U;
            break;
            case 9:	  CHANGE_STAT[j] += edgeflag ? -(double)t030T : (double)t030T;
            break;
            case 10:   CHANGE_STAT[j] += edgeflag ? -(double)t030C : (double)t030C;
            break;
            case 11:  CHANGE_STAT[j] += edgeflag ? -(double)t201 : (double)t201;
            break;
            case 12:  CHANGE_STAT[j] += edgeflag ? -(double)t120D : (double)t120D;
            break;
            case 13:  CHANGE_STAT[j] += edgeflag ? -(double)t120U : (double)t120U;
            break;
            case 14:  CHANGE_STAT[j] += edgeflag ? -(double)t120C : (double)t120C;
            break;
            case 15:  CHANGE_STAT[j] += edgeflag ? -(double)t210 : (double)t210;
            break;
  	        case 16:  CHANGE_STAT[j] += edgeflag ? -(double)t300 : (double)t300;
            break;
          }
        }
        TOGGLE_IF_MORE_TO_COME(i);
    }
  } else {
    /*  undirected */

    /* *** don't forget tail -> head */    
    FOR_EACH_TOGGLE(i) {
      edgeflag = IS_OUTEDGE(tail = TAIL(i), head = HEAD(i));
      t300 = 0; t201 = 0; t102 = 0; t012 = 0;

      if ((EdgetreeMinimum(nwp->outedges, head) != 0) || 
          (EdgetreeMinimum(nwp->inedges, head) != 0) || 
          (EdgetreeMinimum(nwp->outedges, tail) != 0) ||
          (EdgetreeMinimum(nwp->inedges, tail) != 0)) {      

            /* ****** loop through node3 ****** */
            for (node3=1; node3 <= N_NODES; node3++) { 
              if (node3 != tail && node3 != head) {
                a = (EdgetreeSearch(MIN(node3,head), MAX(node3,head), nwp->outedges) != 0);
                b = (EdgetreeSearch(MIN(node3,tail), MAX(node3,tail), nwp->outedges) != 0);
                edgecount = (a + b);
                
                switch(edgecount) {  
                  case 0: {   /* 012 */
                    ++t102;
                    --t012;
                  }
                  break;

                  case 1: {  /* 021C, 021U, 021D, 102 */
                    ++t201;
                    --t102;
                  }
                  break;
          
                  case 2: {  /* 030C, 030T, 111U, 111D */
                    ++t300;
                    --t201;
                  }
                  break;
            
                }
              }

            }    /* ******  move to next node3 ******** */
          } else { 
            t102 = t102 + (N_NODES - 2);  
          }
          
          for(j = 0; j < N_CHANGE_STATS; j++) {
            triadtype = (Vertex)INPUT_PARAM[j]; 
            
            switch(triadtype) { /* SEARCH_ON_THIS_TO_TRACK_DOWN_TRIADCENSUS_CHANGE
                                  to undo triadcensus change, change - to plus in 
                                  next line: */
              case 1:  t003 = -(t102+t201+t300);
              CHANGE_STAT[j] += edgeflag ? -(double)t003 : (double)t003;
              break;
              case 2:  CHANGE_STAT[j] += edgeflag ? -(double)t102 : (double)t102;
              break;
              case 3:  CHANGE_STAT[j] += edgeflag ? -(double)t201 : (double)t201;
              break;
              case 4:  CHANGE_STAT[j] += edgeflag ? -(double)t300 : (double)t300;
              break;
            }
          }
          TOGGLE_IF_MORE_TO_COME(i);
    } /* i loop */
  } 
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_triangle
*****************/
D_CHANGESTAT_FN(d_triangle) { 
  Edge e;
  Vertex tail, head, change, node3;
  int i, j;
  double tailattr, edgemult;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    tail = TAIL(i);
    head = HEAD(i);
    edgemult = IS_OUTEDGE(tail, head) ? -1.0 : 1.0;
    change = 0;
    if(N_INPUT_PARAMS>0){ /* match on attributes */
      tailattr = INPUT_ATTRIB[tail-1];
      if(tailattr == INPUT_ATTRIB[head-1]){
        STEP_THROUGH_OUTEDGES(head, e, node3) { /* step through outedges of head */
          if(tailattr == INPUT_ATTRIB[node3-1]){
            if (DIRECTED) change += IS_OUTEDGE(node3, tail) + IS_INEDGE(node3, tail);
            else change += IS_UNDIRECTED_EDGE(node3,tail);
          }
        }
        STEP_THROUGH_INEDGES(head, e, node3) { /* step through inedges of head */
          if(tailattr == INPUT_ATTRIB[node3-1]){
            if (DIRECTED) change += IS_OUTEDGE(node3, tail) + IS_INEDGE(node3, tail);
            else change += IS_UNDIRECTED_EDGE(node3,tail);
          }
        }
        if(N_CHANGE_STATS>1){ /* diff = TRUE */
          for (j=0; j<N_CHANGE_STATS; j++){
            if (tailattr == INPUT_PARAM[j])
              CHANGE_STAT[j] += edgemult * change;
          }
        }else{ /* diff = FALSE */
          CHANGE_STAT[0] += edgemult * change;
        }
      }
    }else{ /* no attribute matching */
      STEP_THROUGH_OUTEDGES(head, e, node3) { /* step through outedges of head */
        if (DIRECTED) change += IS_OUTEDGE(node3, tail) + IS_INEDGE(node3, tail);
	      else change += IS_UNDIRECTED_EDGE(node3,tail);
      }
      STEP_THROUGH_INEDGES(head, e, node3) { /* step through inedges of head */
        if (DIRECTED) change += IS_OUTEDGE(node3, tail) + IS_INEDGE(node3, tail);
	      else change += IS_UNDIRECTED_EDGE(node3,tail);
      }
      CHANGE_STAT[0] += edgemult * change;
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_tripercent
*****************/
D_CHANGESTAT_FN(d_tripercent) {
  Edge e, e2;
  Vertex tail, head, node1, node2, node3;
  int edgeflag, i, j;
  Edge triwith, triwithout;
  Edge degreewith, degreewithout, twostarwith, twostarwithout;
  int ninputs = N_INPUT_PARAMS - N_NODES;
  int MatchingOnAttribute = (ninputs>0);
  double *attr=INPUT_PARAM, ratiowith, ratiowithout;
  
  if (MatchingOnAttribute) 
    attr = INPUT_PARAM + (ninputs-1); /* ptr to vertex attributes */
 
  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    edgeflag = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i));
    if (!edgeflag) TOGGLE(tail, head); /* turn on the edge if it's missing */
    for (j=0; j < MAX(1, ninputs); j++) {
      /* Count triangles with and without proposed edge */
      /* Simultaneously, find degree (use matching if necessary) with and without */ 
      triwith = triwithout = twostarwith = twostarwithout = 0;
      for (node1 = 1; node1 <= N_NODES; node1++) {
        degreewith = degreewithout = 0;
        if (ninputs < 2 || EQUAL(attr[node1],INPUT_PARAM[j])) {
          STEP_THROUGH_OUTEDGES(node1, e, node2) {
            /* inside this loop, node1 < node2 always */
            if (!MatchingOnAttribute || EQUAL(attr[node1],attr[node2])) {
              /* increment degree counter */
              ++degreewith;
              if (node1!=tail || node2!=head) ++degreewithout;
              STEP_THROUGH_OUTEDGES(node2, e2, node3) {
                /* inside this loop, node1 < node2 < node3 always */
                if (!MatchingOnAttribute || EQUAL(attr[node2],attr[node3])) {
                  if (IS_OUTEDGE(node1, node3)) {
                    ++triwith;
                    if ((tail!=node1||head!=node2)&&(tail!=node2||head!=node3)&&(tail!=node1||head!=node3))
                      ++triwithout;
                  }
                }
              }
            }
          }
          STEP_THROUGH_INEDGES(node1, e, node2) {
            /* inside this loop, node2 < node1 always. */
            /* We only do this to find correct degree for node1; */
            /* for triangles, node1 < node2 <node3 as above suffices */
            if (!MatchingOnAttribute || EQUAL(attr[node1],attr[node2])) {
              /* increment degree counter */
              ++degreewith;
              if (node2!=tail || node1!=head) ++degreewithout;
            }
          }
        }
        /* Now calculate two-star counts (including triangles for now) */
        twostarwith += (degreewith * (degreewith-1))/2; /* same as d-choose-2*/
        twostarwithout += (degreewithout * (degreewithout-1))/2; 
      }
      /* Correct twostar counts for number of triangles */
      twostarwith -= 3*triwith;
      twostarwithout -= 3*triwithout;
      ratiowith = triwith == 0 ? 0.0 : 
                            ((float)triwith)/(triwith + twostarwith);
      ratiowithout = triwith == 0 ? 0.0 :
                            ((float)triwithout)/(triwithout + twostarwithout);
      CHANGE_STAT[j] += (ratiowith-ratiowithout)*(edgeflag? -100.0 : 100.0);
    }
    if (!edgeflag) TOGGLE(tail, head); /* reset dyad to original state */
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_ttriple
*****************/
D_CHANGESTAT_FN(d_ttriple) { 
  Edge e;
  Vertex tail, head, change, node3;
  int i, j;
  double tailattr, edgemult;
  
  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    tail = TAIL(i);
    head = HEAD(i);
    edgemult = IS_OUTEDGE(tail, head) ? -1.0 : 1.0;
    change = 0;
    if(N_INPUT_PARAMS > 0){ /* match on attributes */
      tailattr = INPUT_ATTRIB[tail-1];
      if(tailattr == INPUT_ATTRIB[head-1]) {
        STEP_THROUGH_OUTEDGES(head, e, node3) { /* step through outedges of head */
          if(tailattr == INPUT_ATTRIB[node3-1])
            change += IS_INEDGE(node3, tail);
        }
        STEP_THROUGH_INEDGES(head, e, node3) { /* step through inedges of head */
          if(tailattr == INPUT_ATTRIB[node3-1])
            change += IS_OUTEDGE(node3, tail) + IS_INEDGE(node3, tail);
        }
        if(N_CHANGE_STATS > 1) { /* diff = TRUE; matches must be tabled */
          for (j=0; j<N_CHANGE_STATS; j++){
            if (tailattr == INPUT_PARAM[j])
              CHANGE_STAT[j] += edgemult * change;
          }
        } else { /* diff = FALSE; all matches equivalent */
              CHANGE_STAT[0] += edgemult * change;          
        }
      }
    }else{ /* no attribute matching */
      STEP_THROUGH_OUTEDGES(head, e, node3) { /* step through outedges of head */
        change += IS_INEDGE(node3, tail);
      }
      STEP_THROUGH_INEDGES(head, e, node3) {  /* step through inedges of head */
        change += IS_OUTEDGE(node3, tail) + IS_INEDGE(node3, tail);
      }
      CHANGE_STAT[0] += edgemult * change;
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}


