/*  File src/changestats.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2020 Statnet Commons
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
  double change, absdiff, tailval, headval;
  Vertex tail, head;
  int i, j;
  
  ZERO_ALL_CHANGESTATS(i);

  /* *** don't forget tail -> head */
  FOR_EACH_TOGGLE(i) {
    change = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1.0 : 1.0;
    tailval = INPUT_ATTRIB[tail-1];
    headval = INPUT_ATTRIB[head-1];
    absdiff = fabs(tailval - headval);
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
    echange=IS_OUTEDGE(b1=TAIL(i), HEAD(i)) ? -1:+1;
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
 changestat: d_b1nodematch
*****************/
D_CHANGESTAT_FN(d_b1nodematch) {
  
  Vertex h, t, node3, node4, ninputs;
  int i, edgeflag, count, exponenttype, matchval, b2attrsize, attrval1, attrval2, diffstatus;
  /* int j, numofstats; */
  Edge e, e2;
  double beta, alpha, change=0.0, exponent; 
  const int BetaType=1, AlphaType=2;
  
  b2attrsize = INPUT_PARAM[0];          

  if(b2attrsize > 0){                   
    ninputs = N_INPUT_PARAMS - N_NODES - b2attrsize;/*have 2 sets of node attributes and b2attrvals */  
  }                                      
  else{                                  
    ninputs = N_INPUT_PARAMS - BIPARTITE; 
  }                                      

  diffstatus = !(ninputs == 3); /* 1 if Diff = T and 0 if Diff = F */
  /* numofstats = diffstatus ? (b2attrsize == 0 ? (ninputs - 3): (ninputs - 3) * b2attrsize) : (b2attrsize == 0 ? 1 : b2attrsize); */
  
  exponent = beta = INPUT_PARAM[1]; /* exponent on nodematch count */  
  exponenttype = BetaType;
  alpha = INPUT_PARAM[2];               
  
  if (beta >= 1.0 && alpha < 1.0) {  
    exponent = alpha;
    exponenttype = AlphaType;
  }
  //  Rprintf("N_INPUT_PARAMS = %d, N_NODES=%d\n", N_INPUT_PARAMS, N_NODES);
  //  Rprintf("ninputs = %d, beta=%f, alpha=%f, exponenttype=%d, exponent=%f\n",
  //  ninputs, beta, alpha, exponenttype, exponent);
  
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    t = TAIL(i);
    h = HEAD(i);
    edgeflag = IS_OUTEDGE(t, h);
    matchval = INPUT_PARAM[t + ninputs - 1]; 
    
    /* Now count the neighbors of h whose attribute value equals matchval */
    /* All neighbors of h are inedges because this is a bipartite network */
    count = 0;
    change = 0.0;

    if(b2attrsize == 0){ 
    
      STEP_THROUGH_INEDGES(h, e, node3) {
	    if (INPUT_PARAM[node3 + ninputs - 1] == matchval && t != node3) { /* match! */ 
	        ++count;

	  // Rprintf("Matching twostar found! %d and %d connect to %d\n==================\n", t, node3, h);
	        if (exponenttype == AlphaType) {
	    
	    /* calculate alpha change stat instead of beta change stat. */
	    /* Look for number of two-paths connecting t and node3, not via h */
	        count = 0;
	    
	        STEP_THROUGH_OUTEDGES(t, e2, node4) {
	      // Rprintf("node3=%d, node4=%d, alpha=%f\n", node3,node4,alpha);
		        if (node4 != h) {              /* RPB */
		            count += IS_OUTEDGE(node3, node4); /* add 1 if node4 connects node3 with t */
		        }
	        }
	    
	        /* if count==0, then the statistic is always (plus or minus) 1 */
	        // Rprintf("count is %d\n", count);
	        change += (count==0 ? 1 : pow(count+1, exponent) - pow(count, exponent));
	        }
	      }
        } 
    
      /* If count==0 then the statistic cannot change; it is the same with or */
      /* without the proposed toggle */
    
      if (exponenttype == BetaType && count>0) {
	    /* Now raise count and count+1 to beta, find the difference */
	    change = 0.5*(count+1)*pow(count, exponent);
	    change -= 0.5*count*(exponent==0.0? (count==1? 0.0 : 1.0) : pow((count-1), exponent));
      }

      if (diffstatus) { /* diff=T */
  	    if(ninputs==4) /* keep != NULL*/
	    CHANGE_STAT[0] += edgeflag ? -change : change;
  	    else
	    CHANGE_STAT[matchval-1] += edgeflag ? -change : change;  

      } else { /* diff=F */
	    CHANGE_STAT[0] += edgeflag ? -change : change;
      }

    } else {  
      
      attrval1 = INPUT_PARAM[h + ninputs + b2attrsize - 1];  
 
      STEP_THROUGH_INEDGES(h, e, node3) {
	
	if (INPUT_PARAM[node3 + ninputs - 1] == matchval && t != node3) { /* match! */ 
	 
	  ++count;   

	  // Rprintf("Matching twostar found! %d and %d connect to %d\n==================\n", t, node3, h);
	  if (exponenttype == AlphaType) {
	    /* calculate alpha change stat instead of beta change stat. */
	    /* Look for number of two-paths connecting t and node3, not via h */
	    
	    count = 0;      
	
	    STEP_THROUGH_OUTEDGES(t, e2, node4) {
	      // Rprintf("node3=%d, node4=%d, alpha=%f\n", node3,node4,alpha);
	      if (node4 != h) { 
		    attrval2 = INPUT_PARAM[node4 + ninputs + b2attrsize - 1];  
		    if(attrval2 == attrval1) count += IS_OUTEDGE(node3, node4); 
	      }
	    }
	    /* if count==0, then the statistic is always (plus or minus) 1 */
	    // Rprintf("count is %d\n", count);
	    /* setting the change stat for each parameter */
   	      change += (count== 0 ? 1 : pow(count+1, exponent) - pow(count, exponent)); 	    
        }
	  }
    } 
    /* If count==0 then the statistic cannot change; it is the same with or */
    /* without the proposed toggle */
      if (exponenttype == BetaType && count > 0) {       

      /* Now raise count and count+1 to beta, find the difference */
	   change  = 0.5*(count+1)*pow(count, beta);
	   change -= 0.5*count*(beta==0.0? (count==1? 0.0 : 1.0) : pow((count-1), beta));	
      }       
     
      if(diffstatus){  	 
	      CHANGE_STAT[b2attrsize*(matchval-1) + attrval1 - 1] += edgeflag ? -change : change;
      } else{        
          CHANGE_STAT[attrval1 - 1] += edgeflag ? -change : change;
      }

    } 
    TOGGLE_IF_MORE_TO_COME(i);
  }

  UNDO_PREVIOUS_TOGGLES(i);
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

/*****************
 changestat: d_b2cov
*****************/
D_CHANGESTAT_FN(d_b2cov) { 
  Vertex tail, head;
  int i, edgeflag;
  unsigned int oshift = N_INPUT_PARAMS / N_CHANGE_STATS;

  /* *** don't forget tail -> head */    
  Vertex nb1 = BIPARTITE;
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) 
    {
      edgeflag=IS_OUTEDGE(tail = TAIL(i), head = HEAD(i));
      for(unsigned int j=0, o=0; j<N_CHANGE_STATS; j++, o+=oshift){
	double sum = INPUT_ATTRIB[head-nb1+o-1];
	CHANGE_STAT[j] += edgeflag ? -sum : sum;
      }
      TOGGLE_IF_MORE_TO_COME(i);
    }
  UNDO_PREVIOUS_TOGGLES(i);
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
    echange=IS_OUTEDGE(TAIL(i), b2=HEAD(i)) ? -1:+1;
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
  double s;
  Vertex head;
  int i;
  
  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    head = HEAD(i);
    s = IS_OUTEDGE(TAIL(i), head) ? -1.0 : 1.0;
    int headpos = INPUT_ATTRIB[head-1-BIPARTITE];
    if (headpos!=-1) CHANGE_STAT[headpos] += s;
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_b2nodematch
*****************/
D_CHANGESTAT_FN(d_b2nodematch) {
 
  Vertex h, t, node3, node4, ninputs;
  int i, edgeflag, count, exponenttype, matchval, b1attrsize, attrval1, attrval2, diffstatus;
  /* int j, ind, numofstats; */
  Edge e, e2;
  double beta, alpha, change=0.0, exponent;
  const int BetaType=1, AlphaType=2;

  b1attrsize = INPUT_PARAM[0];             

  if(b1attrsize > 0){                   
    ninputs = N_INPUT_PARAMS - N_NODES - b1attrsize;/*have 2 sets of node attributes and b2attrvals */  
  }                                      
  else{                                  
    ninputs = N_INPUT_PARAMS - N_NODES + BIPARTITE;  
  }

  diffstatus = !(ninputs == 3); /* 1 if Diff = T and o if Diff = F - RPB */
  /* numofstats = diffstatus ? (b1attrsize == 0 ? (ninputs - 3): (ninputs - 3) * b1attrsize) : (b1attrsize == 0 ? 1 : b1attrsize); */

  exponent = beta = INPUT_PARAM[1]; /* exponent on nodematch count */
  exponenttype = BetaType;
  alpha = INPUT_PARAM[2];
  if (beta >= 1.0 && alpha < 1.0) {  
    exponent = alpha;
    exponenttype = AlphaType;
  }
  //  Rprintf("N_INPUT_PARAMS = %d, N_NODES=%d\n", N_INPUT_PARAMS, N_NODES);
  //  Rprintf("ninputs = %d, beta=%f, alpha=%f, exponenttype=%d, exponent=%f\n",
  //  ninputs, beta, alpha, exponenttype, exponent);
  
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    t = TAIL(i);
    h = HEAD(i);
    edgeflag = IS_OUTEDGE(t, h);
    matchval = INPUT_PARAM[h + ninputs - BIPARTITE - 1];
    /* Now count the neighbors of t whose attribute value equals matchval */
    /* All neighbors of t are outedges because this is a bipartite network */
    count=0;
    change = 0.0;
       
    /* RPB */
    /*  double CHANGE[b1attrsize]; */

      
  if(b1attrsize == 0){

    STEP_THROUGH_OUTEDGES(t, e, node3) {
      if (INPUT_PARAM[node3 + ninputs - BIPARTITE - 1] == matchval && h != node3) { /* match! */
        ++count;
        
	// Rprintf("Matching twostar found! %d and %d connect to %d\n==================\n", t, node3, h);
        if (exponenttype == AlphaType) {
          /* calculate alpha change stat instead of beta change stat. */
          /* Look for number of two-paths connecting h and node3 */
          count = 0;
         
      STEP_THROUGH_INEDGES(h, e2, node4) {
            // Rprintf("node3=%d, node4=%d, alpha=%f\n", node3,node4,alpha);
            if (node4 != t) {
              count += IS_OUTEDGE(node4, node3); /* add 1 if node4 connects node3 with h */
            }
          }
          /* if count==0, then the statistic is always 1 */
          // Rprintf("count is %d\n", count);
          change += (count==0 ? 1 : pow(count+1, exponent) - pow(count, exponent));
        }
      }
    }
    /* If count==0 then the statistic cannot change; it is the same with or */
    /* without the proposed toggle */
    
    if (exponenttype == BetaType && count>0) {
      /* Now raise count and count+1 to beta, find the difference */
      change = 0.5*(count+1)*pow(count, beta);
      change -= 0.5*count*(beta==0.0? (count==1? 0.0 : 1.0) : pow((count-1), beta));
    }
    

    if (diffstatus) { /* diff=T */
	    if(ninputs==4) /* keep != NULL*/
	    CHANGE_STAT[0] += edgeflag ? -change : change;
	    else
	    CHANGE_STAT[matchval-1] += edgeflag ? -change : change;

    } else { /* diff=F */
	    CHANGE_STAT[0] += edgeflag ? -change : change;
    }
  } else {

 attrval1 = INPUT_PARAM[t + ninputs + N_NODES + b1attrsize - BIPARTITE - 1];  
 
      STEP_THROUGH_OUTEDGES(t, e, node3) {
	
	if (INPUT_PARAM[node3 + ninputs - BIPARTITE - 1] == matchval && h != node3) { /* match! */ 
	 
	  ++count;   

	  // Rprintf("Matching twostar found! %d and %d connect to %d\n==================\n", t, node3, h);
	  if (exponenttype == AlphaType) {
	    /* calculate alpha change stat instead of beta change stat. */
	    /* Look for number of two-paths connecting t and node3, not via h */
	    
	    count = 0;      
	
	    STEP_THROUGH_INEDGES(h, e2, node4) {
	      // Rprintf("node3=%d, node4=%d, alpha=%f\n", node3,node4,alpha);
	      if (node4 != t) { 
		    attrval2 = INPUT_PARAM[node4 + ninputs + N_NODES + b1attrsize - BIPARTITE - 1];  
		    if(attrval2 == attrval1) count += IS_OUTEDGE(node4, node3); 
	      }
	    }
	    /* if count==0, then the statistic is always (plus or minus) 1 */
	    // Rprintf("count is %d\n", count);
	    /* setting the change stat for each parameter */
   	      change += (count== 0 ? 1 : pow(count+1, exponent) - pow(count, exponent)); 	    
        }
	  }
    } 
    /* If count==0 then the statistic cannot change; it is the same with or */
    /* without the proposed toggle */
      if (exponenttype == BetaType && count > 0) {     

      /* Now raise count and count+1 to beta, find the difference */
	   change  = 0.5*(count+1)*pow(count, beta);
	   change -= 0.5*count*(beta==0.0? (count==1? 0.0 : 1.0) : pow((count-1), beta));	
      }      
     
      if(diffstatus){  	 
	      CHANGE_STAT[b1attrsize*(matchval-1) + attrval1 - 1] += edgeflag ? -change : change;
      } else{        
          CHANGE_STAT[attrval1 - 1] += edgeflag ? -change : change;
      }

     }
    TOGGLE_IF_MORE_TO_COME(i);
  }

  UNDO_PREVIOUS_TOGGLES(i);
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
    STEP_THROUGH_OUTEDGES(head, e, k) /* step through outedges of head */
      {
	if (IS_UNDIRECTED_EDGE(k,tail))
	  ++change;
      }
  }
  
  if(incount){
    STEP_THROUGH_INEDGES(head, e, k) /* step through inedges of head */
      {
	if (IS_UNDIRECTED_EDGE(k,tail))
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
  int i,j,k,semi;
  Vertex tail, head;
  long int maxlen;
  double *countv,emult;
  
  /*Perform initial setup*/
  semi=(int)(INPUT_PARAM[0]);             /*Are we using semicycles?*/
  maxlen=(long int)(INPUT_PARAM[1]);      /*Get max cycle length*/
  countv=Calloc(maxlen-1, double); /*Cycle count vector*/

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    for(j=0;j<maxlen-1;j++)  /*Clear out the count vector*/
      countv[j]=0.0;
    tail = TAIL(i);
    head = HEAD(i);
    /*In semi-cycle case, this toggle can't matter if there is a*/
    /*head->tail edge in the graph; not counting saves much time.*/
    if(!(semi&&(IS_OUTEDGE(head,tail)))){
      /*Count the cycles associated with this edge*/
      edgewise_cycle_census(nwp,tail,head,countv,maxlen,semi);

      /*Make the change, as needed*/
      if((!DIRECTED)&&(tail>head))
        emult = IS_OUTEDGE(head, tail) ? -1.0 : 1.0;
      else
        emult = IS_OUTEDGE(tail, head) ? -1.0 : 1.0;
      k=0;
      for(j=0;j<maxlen-1;j++)
        if(INPUT_PARAM[2+j]>0.0)
          CHANGE_STAT[k++]+=emult*countv[j];
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
  Free(countv);
}

/*****************
 edgewise_path_recurse:  Called by d_cycle
*****************/
void edgewise_path_recurse(Network *nwp, Vertex dest, Vertex curnode, 
     Vertex *visited, long int curlen, double *countv, long int maxlen, int semi) {
  Vertex i,v;
  Edge e;
  int rflag;
  
  /*If we've found a path to the destination, increment the census vector*/ 
  if(DIRECTED){  /*Use outedges, or both if counting semi-paths*/
    if(!semi)
      countv[curlen] += IS_OUTEDGE(curnode, dest);
    else
      countv[curlen] += (IS_OUTEDGE(curnode, dest) || IS_INEDGE(curnode, dest));
  }else{   /*For undirected graphs, edges go from low to high*/
    if(curnode<dest)
      countv[curlen] += IS_OUTEDGE(curnode, dest);
    else
      countv[curlen] += IS_INEDGE(curnode, dest);
  }
  
  /*If possible, keep searching for novel paths*/
  if(curlen<maxlen-2){
    visited[curlen+1]=curnode; /*Add current node to visited list*/

    /*Recurse on all unvisited neighbors of curnode*/
    STEP_THROUGH_OUTEDGES(curnode,e,v){
      rflag=1;
      for(i=0;(i<=curlen)&&(rflag);i++)  /*Check earlier nodes in path*/
        rflag=(v!=visited[i]);
      if(rflag)
        edgewise_path_recurse(nwp,dest,v,visited,curlen+1,countv,maxlen, semi);
    }
    if(semi||(!DIRECTED)){ /*If semi or !directed, need in-neighbors too*/
      STEP_THROUGH_INEDGES(curnode,e,v){
        rflag=((!DIRECTED)||(!(IS_OUTEDGE(curnode,v))));
        for(i=0;(i<=curlen)&&(rflag);i++)  /*Check earlier nodes in path*/
          rflag=(v!=visited[i]);
        if(rflag)
          edgewise_path_recurse(nwp,dest,v,visited,curlen+1,countv,maxlen, semi);
      }
    }
  }
}

/*****************
 edgewise_cycle_census:  Called by d_cycle
*****************/
void edgewise_cycle_census(Network *nwp, Vertex tail, Vertex head, 
                           double *countv, long int maxlen, int semi) {
  /* *** don't forget tail -> head */    
  long int n;
  Vertex *visited,v;
  Edge e;

  /*Set things up*/
  n=N_NODES;

  /*First, check for a 2-cycle (but only if directed and !semi)*/
  if(DIRECTED && (!semi) && IS_OUTEDGE(head,tail))
    countv[0]++;
  if(n==2)
    return;                 /*Failsafe for graphs of order 2*/
  
  /*Perform the recursive path count*/
  visited=Calloc(maxlen,Vertex); /*Initialize the list of visited nodes*/
  visited[0]=tail;
  visited[1]=head;
  
  /*Recurse on each neighbor of head*/
  STEP_THROUGH_OUTEDGES(head,e,v){
    if(v!=tail)
      edgewise_path_recurse(nwp,tail,v,visited,1,countv,maxlen,semi);
  }
  if(semi||(!DIRECTED)){ /*If semi or !directed, need in-neighbors too*/
    STEP_THROUGH_INEDGES(head,e,v){
      if((v!=tail)&&((!DIRECTED)||(!(IS_OUTEDGE(head,v)))))
        edgewise_path_recurse(nwp,tail,v,visited,1,countv,maxlen, semi);
    }
  }
  Free(visited);  /*Free the visited node list*/
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
    echange=IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1:+1;
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
    echange=IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1:+1;
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
  Dyad ndyads = N_DYADS;
  
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
 changestat: d_diff
*****************/
D_CHANGESTAT_FN(d_diff) { 
  double p = INPUT_PARAM[0], *x = INPUT_PARAM+2;
  int mul = INPUT_PARAM[1], sign_code = INPUT_PARAM[2];
  Vertex tail, head;
  int i;

  /* *** don't forget tail -> head */
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    tail = TAIL(i); 
    head = HEAD(i);
    double change = (x[tail] - x[head])*mul;
    switch(sign_code){
    case 1: // identity
      break;
    case 2: // abs
      change = fabs(change);
      break;
    case 3: // positive only
      change = change<0 ? 0 : change;
      break;
    case 4: // negative only
      change = change>0 ? 0 : change;
      break;
    default:
      error("Invalid sign action code passed to d_diff.");
      break;
    }

    if(p==0.0){ // Special case: take the sign of the difference instead.
      change = sign(change);
    }else if(p!=1.0){
      change = pow(change, p);
    }
    
    CHANGE_STAT[0] += IS_OUTEDGE(tail,head) ? -change : change;
    TOGGLE_IF_MORE_TO_COME(i); /* Needed in case of multiple toggles */
  }
  UNDO_PREVIOUS_TOGGLES(i); /* Needed on exit in case of multiple toggles */
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
    echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1:+1;
    /* step through outedges of head */
    STEP_THROUGH_OUTEDGES(head, e, u){
      if (u != tail){
        L2tu=0;
        /* step through outedges of u */
        STEP_THROUGH_OUTEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
        }
        /* step through inedges of u */
        STEP_THROUGH_INEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)INPUT_PARAM[j];
          CHANGE_STAT[j] += ((L2tu + echange == deg)
          - (L2tu == deg));
        }
      }
    }
    /* step through inedges of head */
    STEP_THROUGH_INEDGES(head, e, u){
      if (u != tail){
        L2tu=0;
        /* step through outedges of u */
        STEP_THROUGH_OUTEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
        }
        /* step through inedges of u */
        STEP_THROUGH_INEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)INPUT_PARAM[j];
          CHANGE_STAT[j] += ((L2tu + echange == deg)
          - (L2tu == deg));
        }
      }
    }
    /* step through outedges of tail */
    STEP_THROUGH_OUTEDGES(tail, e, u){
      if (u != head){
        L2uh=0;
        /* step through outedges of u */
        STEP_THROUGH_OUTEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
        }
        /* step through inedges of u */
        STEP_THROUGH_INEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)INPUT_PARAM[j];
          CHANGE_STAT[j] += ((L2uh + echange == deg)
          - (L2uh == deg));
        }
      }
    }
    /* step through inedges of tail */
    STEP_THROUGH_INEDGES(tail, e, u){
      if (u != head){
        L2uh=0;
        /* step through outedges of u */
        STEP_THROUGH_OUTEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
        }
        /* step through inedges of u */
        STEP_THROUGH_INEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
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
      refedgeflag = (IS_OUTEDGE(head, tail));
      
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
D_CHANGESTAT_FN(d_edges){
  int edgeflag, i;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS();
  FOR_EACH_TOGGLE(i){
    edgeflag = IS_OUTEDGE(TAIL(i), HEAD(i));
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
    echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1:+1;
    /* step through outedges of head */
    STEP_THROUGH_OUTEDGES(head, e, u) {
      if (IS_UNDIRECTED_EDGE(u, tail)){
        L2th++;
        L2tu=0;
        L2uh=0;
        /* step through outedges of u */
        STEP_THROUGH_OUTEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
          if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
        }
        /* step through inedges of u */
        STEP_THROUGH_INEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
          if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
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
    STEP_THROUGH_INEDGES(head, e, u){
      if (IS_UNDIRECTED_EDGE(u, tail)){
        L2th++;
        L2tu=0;
        L2uh=0;
        /* step through outedges of u */
        STEP_THROUGH_OUTEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
          if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
        }
        /* step through inedges of u */
        STEP_THROUGH_INEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
          if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
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
    echange=IS_OUTEDGE(b1=TAIL(i), HEAD(i)) ? -1 : +1;
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
    echange=IS_OUTEDGE(b1=TAIL(i), HEAD(i)) ? -1 : +1;
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
    echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1:+1;
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
    echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1:+1;
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
    ochange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 0;
    echange = 2*ochange + 1;
    /* step through outedges of head */
    STEP_THROUGH_OUTEDGES(head, e, u){
      if (u != tail){
        L2tu=ochange;
        /* step through outedges of u */
        STEP_THROUGH_OUTEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
        }
        /* step through inedges of u */
        STEP_THROUGH_INEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
        }
        cumchange += pow(oneexpa,(double)L2tu);
      }
    }
    /* step through inedges of head */
    STEP_THROUGH_INEDGES(head, e, u){
      if (u != tail){
        L2tu=ochange;
        /* step through outedges of u */
        STEP_THROUGH_OUTEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
        }
        /* step through inedges of u */
        STEP_THROUGH_INEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
        }
        cumchange += pow(oneexpa,(double)L2tu);
      }
    }
    
    /* step through outedges of tail  */
    STEP_THROUGH_OUTEDGES(tail, e, u){
      if (u != head){
        L2uh=ochange;
        /* step through outedges of u */
        STEP_THROUGH_OUTEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
        }
        /* step through inedges of u */
        STEP_THROUGH_INEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
        }
        cumchange += pow(oneexpa,(double)L2uh);
      }
    }
    /* step through inedges of tail */
    STEP_THROUGH_INEDGES(tail, e, u){
      if (u != head){
        L2uh=ochange;
        /* step through outedges of u */
        STEP_THROUGH_OUTEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
        }
        /* step through inedges of u */
        STEP_THROUGH_INEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
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
    echange=IS_OUTEDGE(TAIL(i), b2=HEAD(i)) ? -1 : +1;
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
    echange=IS_OUTEDGE(TAIL(i), b2=HEAD(i)) ? -1 : +1;
    b2deg = id[b2]+(echange-1)/2;
    b2attr = INPUT_PARAM[b2 - BIPARTITE]; 
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
    ochange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 0;
    echange = 2*ochange + 1;
    /* step through outedges of head  */
    STEP_THROUGH_OUTEDGES(head, e, u){
      if (IS_UNDIRECTED_EDGE(u, tail)){
	L2th++;
	L2tu=ochange;
	L2uh=ochange;
	/* step through outedges of u */
	STEP_THROUGH_OUTEDGES(u, f, v){
	  if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
	  if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
	}
	/* step through inedges of u */
	STEP_THROUGH_INEDGES(u, f, v){
	  if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
	  if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
	}
	cumchange += pow(oneexpa,(double)L2tu) +
	  pow(oneexpa,(double)L2uh) ;
      }
    }
    /* step through inedges of head */
    
    STEP_THROUGH_INEDGES(head, e, u){
      if (IS_UNDIRECTED_EDGE(u, tail)){
	L2th++;
	L2tu=ochange;
	L2uh=ochange;
	/* step through outedges of u */
	STEP_THROUGH_OUTEDGES(u, f, v){
	  if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
	  if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
	}
	/* step through inedges of u */
	STEP_THROUGH_INEDGES(u, f, v){
	  if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
	  if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
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
    ochange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 0;
    echange = 2*ochange + 1;
    /* step through outedges of head */
    STEP_THROUGH_OUTEDGES(head, e, u){
      if (u != tail){
        L2tu=ochange;
        /* step through outedges of u */
        STEP_THROUGH_OUTEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
        }
        /* step through inedges of u */
        STEP_THROUGH_INEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
        }
        cumchange += pow(oneexpa,(double)L2tu);
      }
    }
    /* step through inedges of head */
    STEP_THROUGH_INEDGES(head, e, u){
      if (u != tail){
        L2tu=ochange;
        /* step through outedges of u */
        STEP_THROUGH_OUTEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
        }
        /* step through inedges of u */
        STEP_THROUGH_INEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
        }
        cumchange += pow(oneexpa,(double)L2tu);
      }
    }
    
    /* step through outedges of tail  */
    STEP_THROUGH_OUTEDGES(tail, e, u){
      if (u != head){
        L2uh=ochange;
        /* step through outedges of u */
        STEP_THROUGH_OUTEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
        }
        /* step through inedges of u */
        STEP_THROUGH_INEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
        }
        cumchange += pow(oneexpa,(double)L2uh);
      }
    }
    /* step through inedges of tail */
    STEP_THROUGH_INEDGES(tail, e, u){
      if (u != head){
        L2uh=ochange;
        /* step through outedges of u */
        STEP_THROUGH_OUTEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
        }
        /* step through inedges of u */
        STEP_THROUGH_INEDGES(u, f, v){
          if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
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
    ochange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 0;
    echange = 2*ochange + 1;
    /* step through outedges of head  */
    STEP_THROUGH_OUTEDGES(head, e, u){
      if (IS_UNDIRECTED_EDGE(u, tail)){
	L2th++;
	L2tu=ochange;
	L2uh=ochange;
	/* step through outedges of u */
	STEP_THROUGH_OUTEDGES(u, f, v){
	  if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
	  if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
	}
	/* step through inedges of u */
	STEP_THROUGH_INEDGES(u, f, v){
	  if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
	  if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
	}
	cumchange += pow(oneexpa,(double)L2tu) +
	  pow(oneexpa,(double)L2uh) ;
      }
    }
    /* step through inedges of head */
    
    STEP_THROUGH_INEDGES(head, e, u){
      if (IS_UNDIRECTED_EDGE(u, tail)){
	L2th++;
	L2tu=ochange;
	L2uh=ochange;
	/* step through outedges of u */
	STEP_THROUGH_OUTEDGES(u, f, v){
	  if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
	  if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
	}
	/* step through inedges of u */
	STEP_THROUGH_INEDGES(u, f, v){
	  if(IS_UNDIRECTED_EDGE(v, head)) L2uh++;
	  if(IS_UNDIRECTED_EDGE(v, tail)) L2tu++;
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
    ochange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 0;
    echange = 2*ochange + 1;
    /* step through outedges of head  */
    STEP_THROUGH_OUTEDGES(head, e, u){
      if (IS_OUTEDGE(tail, u)){
	L2tu=ochange;
	/* step through inedges of u */
	STEP_THROUGH_INEDGES(u, f, v){
	  if(IS_OUTEDGE(tail, v)) L2tu++;
	}
	cumchange += pow(oneexpa,(double)L2tu);
      }
    }
    /* step through inedges of head */
    
    STEP_THROUGH_INEDGES(head, e, u){
      if (IS_OUTEDGE(tail, u)){
	L2th++;
      }
      if (IS_OUTEDGE(u, tail)){
	L2uh=ochange;
	/* step through outedges of u */
	STEP_THROUGH_OUTEDGES(u, f, v){
	  if(IS_OUTEDGE(v, head)) L2uh++;
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
    ochange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 0;
    echange = 2*ochange + 1;
    /* step through outedges of head  */
    STEP_THROUGH_OUTEDGES(head, e, u){
      if (IS_OUTEDGE(tail, u)){
	L2tu=ochange;
	/* step through inedges of u */
	STEP_THROUGH_INEDGES(u, f, v){
	  if(IS_OUTEDGE(tail, v)) L2tu++;
	}
	cumchange += pow(oneexpa,(double)L2tu);
      }
    }
    /* step through inedges of head */
    
    STEP_THROUGH_INEDGES(head, e, u){
      if (IS_OUTEDGE(tail, u)){
	L2th++;
      }
      if (IS_OUTEDGE(u, tail)){
	L2uh=ochange;
	/* step through outedges of u */
	STEP_THROUGH_OUTEDGES(u, f, v){
	  if(IS_OUTEDGE(v, head)) L2uh++;
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
    echange=IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1:+1;
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
    int echange = IS_OUTEDGE(TAIL(i), head=HEAD(i)) ? -1 : +1;
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
    echange=IS_OUTEDGE(TAIL(i), head=HEAD(i)) ? -1 : +1;
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
      echange=IS_OUTEDGE(tail, head) ? -1 : +1;
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
        if (!IS_OUTEDGE(tail,node2)){
          change = change + 1.0;
        }
      }
    }
    STEP_THROUGH_INEDGES(head, e, node2) {
      if (node2 != tail){
        if (IS_OUTEDGE(tail, node2)){
          change = change - 1.0;
        }
      }
    }
    STEP_THROUGH_INEDGES(tail, e, node2) {
      if (node2 != head){
        if (!IS_OUTEDGE(node2,head)){
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
changestat: d_isolatededges
*****************/
D_CHANGESTAT_FN(d_isolatededges) { 
  int i, edgeflag;
  Vertex tail, head, neighbor, taild, headd, *id, *od;
  Edge e;
  
  id=IN_DEG;
  od=OUT_DEG;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    // is there an edge tail -> head?
    edgeflag = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i));
    
    taild = od[tail] + id[tail];
    headd = od[head] + id[head];
    
    if(edgeflag) { // we are removing an edge
        
      // if head and tail both have degree one, then
      // we are removing an isolated edge        
      if(taild == 1 && headd == 1)
        CHANGE_STAT[0] -= 1;
          
      // if tail has degree 2 and has a degree one node other than head as a neighbor,
      // then we are making a non-isolated edge into an isolated edge by removing
      // the edge tail -> head
      if(taild == 2) { 
        STEP_THROUGH_OUTEDGES(tail, e, neighbor) {
          if(od[neighbor] + id[neighbor] == 1 && neighbor != head)
            CHANGE_STAT[0] += 1;
        }
        STEP_THROUGH_INEDGES(tail, e, neighbor) {
          if(od[neighbor] + id[neighbor] == 1 && neighbor != head)
            CHANGE_STAT[0] += 1;
        }          
      }
        
      // ditto head
      if(headd == 2) { 
        STEP_THROUGH_OUTEDGES(head, e, neighbor) {
          if(od[neighbor] + id[neighbor] == 1 && neighbor != tail)
            CHANGE_STAT[0] += 1;
        }
        STEP_THROUGH_INEDGES(head, e, neighbor) {
          if(od[neighbor] + id[neighbor] == 1 && neighbor != tail)
            CHANGE_STAT[0] += 1;
        }        
      }
    } else { // we are adding an edge
    
      // if head and tail both have degree zero, then
      // we are adding an isolated edge
      if(taild == 0 && headd == 0)
        CHANGE_STAT[0] += 1;
      
      // if tail has degree 1 and so does its neighbor, then we 
      // are making an isolated edge into a non-isolated edge;
      // note that for undirected graphs, this neighbor cannot
      // be head, as the current toggle is to turn on the edge
      // tail -> head
      if(taild == 1) { 
        STEP_THROUGH_OUTEDGES(tail, e, neighbor) {
          if(od[neighbor] + id[neighbor] == 1)
            CHANGE_STAT[0] -= 1;
        }
        STEP_THROUGH_INEDGES(tail, e, neighbor) {
          if(od[neighbor] + id[neighbor] == 1)
            CHANGE_STAT[0] -= 1;
        }
      }
      
      // ditto head
      if(headd == 1) { 
        STEP_THROUGH_OUTEDGES(head, e, neighbor) {
          if(od[neighbor] + id[neighbor] == 1)
            CHANGE_STAT[0] -= 1;
        }
        STEP_THROUGH_INEDGES(head, e, neighbor) {
          if(od[neighbor] + id[neighbor] == 1)
            CHANGE_STAT[0] -= 1;
        }
      }
    }
    
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
      echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1:+1;
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
        STEP_THROUGH_INEDGES(head, e, node3) {/* step through inedges of head */
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
      backedgeflag = (IS_OUTEDGE(head, tail));

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
    CHANGE_STAT[0] += (IS_OUTEDGE(TAIL(i), HEAD(i)) ? -2 : 2);
    TOGGLE_IF_MORE_TO_COME(i);
  }
  CHANGE_STAT[0]/=N_NODES*(DIRECTED+1); // Effectively, change is 2/n if undirected and 1/n if directed.
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
 changestat: d_mixmat
 General mixing matrix (mm) implementation.
*****************/
D_CHANGESTAT_FN(d_mixmat){
  unsigned int symm = ((int)INPUT_PARAM[0]) & 1;
  unsigned int marg = ((int)INPUT_PARAM[0]) & 2;
  double *tx = INPUT_PARAM;
  double *hx = BIPARTITE? INPUT_PARAM : INPUT_PARAM + N_NODES;
  double *cells = BIPARTITE? INPUT_PARAM + N_NODES + 1: INPUT_PARAM + N_NODES*2 + 1;
  
  int i;
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i){
    Vertex tail=TAIL(i);
    Vertex head=HEAD(i);
    unsigned int edgeflag = IS_OUTEDGE(tail, head);
    unsigned int diag = tx[tail]==tx[head] && hx[tail]==hx[head];
    for(unsigned int j=0; j<N_CHANGE_STATS; j++){
      unsigned int thmatch = tx[tail]==cells[j*2] && hx[head]==cells[j*2+1];
      unsigned int htmatch = tx[head]==cells[j*2] && hx[tail]==cells[j*2+1];
      
      int w = DIRECTED || BIPARTITE? thmatch :
	(symm ? thmatch||htmatch : thmatch+htmatch)*(symm && marg && diag?2:1);
      if(w) CHANGE_STAT[j] += edgeflag ? -w : w;
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
  edgeflagth = (!IS_OUTEDGE(head,tail));
   
  for(node3=1;node3<=N_NODES;node3++){
    if((node3!=tail)&&(node3!=head)){
     sc = edgeflagth + (!IS_OUTEDGE(node3,tail));
     if(sc < 2){
      sc += (!IS_OUTEDGE(tail,node3));
      if(sc < 2){
       sc += (!IS_OUTEDGE(node3,head));
       if(sc < 2){
        sc += (!IS_OUTEDGE(head,node3));
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
  Vertex tail, head;
  int i, edgeflag;
  unsigned int oshift = N_INPUT_PARAMS / N_CHANGE_STATS;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) 
    {
      edgeflag=IS_OUTEDGE(tail = TAIL(i), head = HEAD(i));
      for(unsigned int j=0, o=0; j<N_CHANGE_STATS; j++, o+=oshift){
	double sum = INPUT_ATTRIB[tail+o-1] + INPUT_ATTRIB[head+o-1];
	CHANGE_STAT[j] += edgeflag ? -sum : sum;
      }
      TOGGLE_IF_MORE_TO_COME(i);
    }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_nodefactor
*****************/
D_CHANGESTAT_FN(d_nodefactor) { 
  double s;
  Vertex tail, head;
  int i;
  
  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    tail = TAIL(i);
    head = HEAD(i);
    s = IS_OUTEDGE(tail, head) ? -1.0 : 1.0;
    int tailpos = INPUT_ATTRIB[tail-1];
    int headpos = INPUT_ATTRIB[head-1];
    if (tailpos!=-1) CHANGE_STAT[tailpos] += s;
    if (headpos!=-1) CHANGE_STAT[headpos] += s;
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_nodeicov
*****************/
D_CHANGESTAT_FN(d_nodeicov) { 
  Vertex tail, head;
  int i, edgeflag;
  unsigned int oshift = N_INPUT_PARAMS / N_CHANGE_STATS;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) 
    {
      edgeflag=IS_OUTEDGE(tail = TAIL(i), head = HEAD(i));
      for(unsigned int j=0, o=0; j<N_CHANGE_STATS; j++, o+=oshift){
	double sum = INPUT_ATTRIB[head+o-1];
	CHANGE_STAT[j] += edgeflag ? -sum : sum;
      }
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
  int i;
  
  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    head = HEAD(i);
    s = IS_OUTEDGE(TAIL(i), head) ? -1.0 : 1.0;
    int headpos = INPUT_ATTRIB[head-1];
    if (headpos!=-1) CHANGE_STAT[headpos] += s;
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
  Vertex tail, head;
  int i, edgeflag;
  unsigned int oshift = N_INPUT_PARAMS / N_CHANGE_STATS;

  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) 
    {
      edgeflag=IS_OUTEDGE(tail = TAIL(i), head = HEAD(i));
      for(unsigned int j=0, o=0; j<N_CHANGE_STATS; j++, o+=oshift){
	double sum = INPUT_ATTRIB[tail+o-1];
	CHANGE_STAT[j] += edgeflag ? -sum : sum;
      }
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
  int i;
  
  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    tail = TAIL(i);
    s = IS_OUTEDGE(tail, HEAD(i)) ? -1.0 : 1.0;
    int tailpos = INPUT_ATTRIB[tail-1];
    if (tailpos!=-1) CHANGE_STAT[tailpos] += s;
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
    echange=IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1:+1;
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
    echange=IS_OUTEDGE(tail=TAIL(i), HEAD(i)) ? -1 : +1;
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
      int echange=IS_OUTEDGE(tail, head) ? -1 : +1;
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
        STEP_THROUGH_OUTEDGES(tail, e, node3) { /* step through outedges of head */
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

/********************  changestats:  P    ***********/
/*****************
 changestat: d_pdegcor
*****************/
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

/********************  changestats:  R    ***********/
/*****************
 changestat: d_rdegcor
*****************/
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
   
  if(IS_OUTEDGE(head, tail)){
   change = 0;
   
   STEP_THROUGH_OUTEDGES(head, e, node3) /* step through outedges of head */
   {
     if (node3 != tail
      && IS_OUTEDGE(node3, tail) 
      && IS_OUTEDGE(tail, node3) 
      && IS_OUTEDGE(node3, head) 
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
    ((IS_OUTEDGE(tail, head)) ? -1.0 : 1.0); 
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
    echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1:+1;
    /* step through outedges of head */
    STEP_THROUGH_OUTEDGES(head, e, u) {
      if (IS_OUTEDGE(tail, u)){
        L2tu=0;
        /* step through inedges of u */
        STEP_THROUGH_INEDGES(u, f, v){
          if(IS_OUTEDGE(tail, v)) L2tu++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)INPUT_PARAM[j];
          CHANGE_STAT[j] += ((L2tu + echange == deg) - (L2tu == deg));
        }
      }
    }
    /* step through inedges of head */
    STEP_THROUGH_INEDGES(head, e, u){
      if (IS_OUTEDGE(tail, u)){
        L2th++;
      }
      if (IS_OUTEDGE(u, tail)){
        L2uh=0;
        /* step through outedges of u */
        STEP_THROUGH_OUTEDGES(u, f, v){
          if(IS_OUTEDGE(v, head)) L2uh++;
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
 changestat: d_threetrail
*****************/
D_CHANGESTAT_FN(d_threetrail) { 
  int i, j, k, edgeflag, change, dchange[4];
  Edge e;
  Vertex tail, head, node3;
  /* The four values of dchange represent the four different types of
     directed threetrails oriented so that the middle step is always 
     "right" (R).  In order:   RRR, RRL, LRR, LRL 
                         i.e., >>>  >><  <>>  <><  */


  /* *** don't forget tail -> head */    
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    edgeflag = IS_OUTEDGE(tail = TAIL(i), head = HEAD(i));
    /* Step A: Count threetrails in which tail->head is the middle edge */
    dchange[0] = IN_DEG[tail] * OUT_DEG[head]; /* R then R; may count head->tail->head->tail */
    dchange[1] = IN_DEG[tail] * (IN_DEG[head]-edgeflag); /* R then L */
    dchange[2] = (OUT_DEG[tail]-edgeflag) * OUT_DEG[head]; /* L then R */
    dchange[3] = (OUT_DEG[tail]-edgeflag) * (IN_DEG[head]-edgeflag); /* L then L */
    /* Step B: Count threetrails where tail is one endpoint */
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
    /* Step C: Count threetrails where head is one endpoint */
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
    echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1:+1;
    /* step through outedges of head */
    STEP_THROUGH_OUTEDGES(head, e, u) {
      if (IS_OUTEDGE(tail, u)){
        L2tu=0;
        /* step through inedges of u */
        STEP_THROUGH_INEDGES(u, f, v){
          if(IS_OUTEDGE(tail, v)) L2tu++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)INPUT_PARAM[j];
          CHANGE_STAT[j] -= ((L2tu + echange == deg) - (L2tu == deg));
        }
      }
    }
    /* step through inedges of head */
    STEP_THROUGH_INEDGES(head, e, u){
      if (IS_OUTEDGE(tail, u)){
        L2th++;
      }
      if (IS_OUTEDGE(u, tail)){
        L2uh=0;
        /* step through outedges of u */
        STEP_THROUGH_OUTEDGES(u, f, v){
          if(IS_OUTEDGE(v, head)) L2uh++;
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
    ochange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 0;
    echange = 2*ochange + 1;
    if(N_INPUT_PARAMS>0){ /* match on attributes */
      tailattr = INPUT_ATTRIB[tail-1];
      if(tailattr == INPUT_ATTRIB[head-1]){
       /* step through outedges of head  */
       STEP_THROUGH_OUTEDGES(head, e, u){
         if (IS_OUTEDGE(tail, u) && (tailattr == INPUT_ATTRIB[u-1])){
	   L2tu=ochange;
	   /* step through inedges of u */
	   STEP_THROUGH_INEDGES(u, f, v){
	     if(IS_OUTEDGE(tail, v) && (tailattr == INPUT_ATTRIB[v-1])){
	       L2tu++;
	       if(L2tu>0) {break;}
	     }
	   }
	   cumchange += (L2tu==0);
         }
       }
       /* step through inedges of head */
       
       STEP_THROUGH_INEDGES(head, e, u){
         if (IS_OUTEDGE(tail, u) && (tailattr == INPUT_ATTRIB[u-1])){
	   L2th++;
         }
         if (IS_OUTEDGE(u, tail) && (tailattr == INPUT_ATTRIB[u-1])){
	   L2uh=ochange;
	   /* step through outedges of u */
	   STEP_THROUGH_OUTEDGES(u, f, v){
	     if(IS_OUTEDGE(v, head) && (tailattr == INPUT_ATTRIB[v-1])){
	       L2uh++;
	       if(L2uh>0) {break;}
	     }
	   }
	   cumchange += (L2uh==0) ;
         }
       }}
      }else{ /* no attributes */
    /* step through outedges of head  */
    STEP_THROUGH_OUTEDGES(head, e, u){
      if (IS_OUTEDGE(tail, u)){
	L2tu=ochange;
	/* step through inedges of u */
	STEP_THROUGH_INEDGES(u, f, v){
	  if(IS_OUTEDGE(tail, v)){
	    L2tu++;
	    if(L2tu>0) {break;}
	  }
	}
	cumchange += (L2tu==0);
      }
    }
    /* step through inedges of head */
    
    STEP_THROUGH_INEDGES(head, e, u){
      if (IS_OUTEDGE(tail, u)){
	L2th++;
      }
      if (IS_OUTEDGE(u, tail)){
	L2uh=ochange;
	/* step through outedges of u */
	STEP_THROUGH_OUTEDGES(u, f, v){
	  if(IS_OUTEDGE(v, head)){
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
    ochange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 0;
    echange = 2*ochange + 1;
    if(N_INPUT_PARAMS>0){ /* match on attributes */
      tailattr = INPUT_ATTRIB[tail-1];
      if(tailattr == INPUT_ATTRIB[head-1]){
       /* step through outedges of head  */
       STEP_THROUGH_OUTEDGES(head, e, u){
         if (IS_INEDGE(tail, u) && (tailattr == INPUT_ATTRIB[u-1])){
	   L2tu=ochange;
	   /* step through inedges of u */
	   STEP_THROUGH_INEDGES(u, f, v){
	     if(IS_OUTEDGE(tail, v) && (tailattr == INPUT_ATTRIB[v-1])){
	       L2tu++;
	       if(L2tu>0) {break;}
	     }
	   }
	   cumchange += (L2tu==0);
         }
       }
       /* step through inedges of head */
       
       STEP_THROUGH_OUTEDGES(head, e, u){
         if (IS_INEDGE(tail, u) && (tailattr == INPUT_ATTRIB[u-1])){
	   L2th++;
         }
         if (IS_OUTEDGE(u, tail) && (tailattr == INPUT_ATTRIB[u-1])){
	   L2uh=ochange;
	   /* step through outedges of u */
	   STEP_THROUGH_OUTEDGES(u, f, v){
	     if(IS_OUTEDGE(v, head) && (tailattr == INPUT_ATTRIB[v-1])){
	       L2uh++;
	       if(L2uh>0) {break;}
	     }
	   }
	   cumchange += (L2uh==0) ;
         }
       }}
      }else{ /* no attributes */
    /* step through outedges of head  */
    STEP_THROUGH_OUTEDGES(head, e, u){
      if (IS_INEDGE(tail, u)){
	L2tu=ochange;
	/* step through inedges of u */
	STEP_THROUGH_INEDGES(u, f, v){
	  if(IS_OUTEDGE(tail, v)){
	    L2tu++;
	    if(L2tu>0) {break;}
	  }
	}
	cumchange += (L2tu==0);
      }
    }
    /* step through outedges of head */
    
    STEP_THROUGH_OUTEDGES(head, e, u){
      if (IS_INEDGE(tail, u)){
	L2th++;
      }
      if (IS_OUTEDGE(u, tail)){
	L2uh=ochange;
	/* step through outedges of u */
	STEP_THROUGH_OUTEDGES(u, f, v){
	  if(IS_OUTEDGE(v, head)){
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
      t120C = 0; t120U = 0; t120D = 0;  t201 = 0;
      t030C = 0; t030T = 0; t111U = 0; t111D = 0;
      t021C = 0; t021U = 0; t021D = 0;  t102 = 0;
       t012 = 0;
      
      if ( (MIN_OUTEDGE(head) != 0) || 
           (MIN_INEDGE(head)  != 0) || 
           (MIN_OUTEDGE(tail) != 0) ||
           (MIN_INEDGE(tail)  != 0) ) {      

          /* ****** loop through node3 ****** */
          for (node3=1; node3 <= N_NODES; node3++) { 
            if (node3 != tail && node3 != head) {
              a = (IS_OUTEDGE(head, tail)); 
              b = (IS_OUTEDGE(head, node3));
              c = (IS_OUTEDGE(node3, head));
              d = (IS_OUTEDGE(node3, tail));
              e = (IS_OUTEDGE(tail, node3));
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
        }else{
          t012 = t012 + (N_NODES - 2);  
	}

        for(j = 0; j < N_CHANGE_STATS; j++) { 
          triadtype = (Vertex)INPUT_PARAM[j]; 
          
          switch(triadtype) { /* SEARCH_ON_THIS_TO_TRACK_DOWN_TRIADCENSUS_CHANGE
                                 to undo triadcensus change, change - to plus in 
                                  next two lines: */
            case 1:   t003 = -(t300+t210+t120C+t120U+t120D+t201+t030C+t030T);
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
            case 7:   CHANGE_STAT[j] += edgeflag ? -(double)t111D : (double)t111D;
            break;
            case 8:   CHANGE_STAT[j] += edgeflag ? -(double)t111U : (double)t111U;
            break;
            case 9:   CHANGE_STAT[j] += edgeflag ? -(double)t030T : (double)t030T;
            break;
            case 10:  CHANGE_STAT[j] += edgeflag ? -(double)t030C : (double)t030C;
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

      if ( (MIN_OUTEDGE(head) != 0) || 
           (MIN_INEDGE(head)  != 0) || 
           (MIN_OUTEDGE(tail) != 0) ||
           (MIN_INEDGE(tail)  != 0) ) {      

            /* ****** loop through node3 ****** */
            for (node3=1; node3 <= N_NODES; node3++) { 
              if (node3 != tail && node3 != head) {
                a = (IS_UNDIRECTED_EDGE(node3, head));
                b = (IS_UNDIRECTED_EDGE(node3, tail));
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


