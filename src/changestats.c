/*  File src/changestats.c in package ergm, part of the Statnet suite of
 *  packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#include "changestats.h"
#include "ergm_storage.h"
#include "ergm_dyad_hashmap.h"
#include "ergm_edgelist.h"

/********************  changestats:  A    ***********/

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
    taildeg = DEG(tail);
    headdeg = DEG(head);
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
C_CHANGESTAT_FN(c_altkstar) {
  double lambda, oneexpl, change;
  Vertex taild, headd=0;

  change = 0.0;
  lambda = INPUT_PARAM[0];
  oneexpl = 1.0-1.0/lambda;

  /* *** don't forget tail -> head */
    taild = DEG(tail) - edgestate;
    headd = DEG(head) - edgestate;
    if(taild!=0){
      change += (edgestate?-1:+1)*(1.0-pow(oneexpl,(double)taild));
    }
    if(headd!=0){
      change += (edgestate?-1:+1)*(1.0-pow(oneexpl,(double)headd));
    }
  CHANGE_STAT[0] = change*lambda;
}

/*****************
 changestat: d_asymmetric
*****************/
C_CHANGESTAT_FN(c_asymmetric) {
  double matchval, change;
  int j, ninputs, noattr;

  ninputs = N_INPUT_PARAMS - N_NODES;
  noattr = (N_INPUT_PARAMS == 0);

  /* *** don't forget tail -> head */
    change = (edgestate==IS_OUTEDGE(head, tail) ? 1.0 : -1.0) ;
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
C_CHANGESTAT_FN(c_b1concurrent) {
  int echange;
  Vertex b1, b1deg;

  /* *** don't forget tail -> head */
    b1 = tail;
    echange = IS_OUTEDGE(b1,head) ? -1 : 1;
    b1deg = OUT_DEG[b1];
    CHANGE_STAT[0] += (b1deg + echange > 1) - (b1deg > 1);
}

/*****************
 changestat: d_b1concurrent_by_attr
*****************/
C_CHANGESTAT_FN(c_b1concurrent_by_attr) {
  /*The inputparams are assumed to be set up as follows:
    The first 2*nstats values are in pairs:  (degree, attrvalue)
    The values following the first 2*nstats values are the nodal attributes. */
  int j, echange, b1attr;
  Vertex b1, b1deg;

  /* *** don't forget tail -> head */
    b1 = tail;
    echange = IS_OUTEDGE(b1,head) ? -1 : 1;
    b1deg = OUT_DEG[b1];
    b1attr = INPUT_PARAM[N_CHANGE_STATS + b1 - 1];
    for(j = 0; j < N_CHANGE_STATS; j++) {
      if (b1attr == INPUT_PARAM[j]) { /* we have attr match */
        CHANGE_STAT[j] += (b1deg + echange > 1) - (b1deg > 1);
      }
    }
}

/*****************
 changestat: d_b1nodematch
*****************/
C_CHANGESTAT_FN(c_b1nodematch) {

  Vertex node3, node4, ninputs;
  int count, exponenttype, matchval, b2attrsize, attrval1, attrval2, diffstatus;
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

    matchval = INPUT_PARAM[tail + ninputs - 1];

    /* Now count the neighbors of head whose attribute value equals matchval */
    /* All neighbors of head are inedges because this is a bipartite network */
    count = 0;
    change = 0.0;

    if(b2attrsize == 0){

      STEP_THROUGH_INEDGES(head, e, node3) {
	    if (INPUT_PARAM[node3 + ninputs - 1] == matchval && tail != node3) { /* match! */
	        ++count;

	  // Rprintf("Matching twostar found! %d and %d connect to %d\n==================\n", tail, node3, head);
	        if (exponenttype == AlphaType) {

	    /* calculate alpha change stat instead of beta change stat. */
	    /* Look for number of two-paths connecting tail and node3, not via head */
	        count = 0;

	        STEP_THROUGH_OUTEDGES(tail, e2, node4) {
	      // Rprintf("node3=%d, node4=%d, alpha=%f\n", node3,node4,alpha);
		        if (node4 != head) {              /* RPB */
		            count += IS_OUTEDGE(node3, node4); /* add 1 if node4 connects node3 with tail */
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
	    CHANGE_STAT[0] += edgestate ? -change : change;
	    else
	    CHANGE_STAT[matchval-1] += edgestate ? -change : change;

      } else { /* diff=F */
	    CHANGE_STAT[0] += edgestate ? -change : change;
      }

    } else {

      attrval1 = INPUT_PARAM[head + ninputs + b2attrsize - 1];

      STEP_THROUGH_INEDGES(head, e, node3) {

	if (INPUT_PARAM[node3 + ninputs - 1] == matchval && tail != node3) { /* match! */

	  ++count;

	  // Rprintf("Matching twostar found! %d and %d connect to %d\n==================\n", tail, node3, head);
	  if (exponenttype == AlphaType) {
	    /* calculate alpha change stat instead of beta change stat. */
	    /* Look for number of two-paths connecting tail and node3, not via head */

	    count = 0;

	    STEP_THROUGH_OUTEDGES(tail, e2, node4) {
	      // Rprintf("node3=%d, node4=%d, alpha=%f\n", node3,node4,alpha);
	      if (node4 != head) {
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
	      CHANGE_STAT[b2attrsize*(matchval-1) + attrval1 - 1] += edgestate ? -change : change;
      } else{
          CHANGE_STAT[attrval1 - 1] += edgestate ? -change : change;
      }

    }
}

/*****************
 changestat: d_b1starmix
*****************/
C_CHANGESTAT_FN(c_b1starmix) {
  double change;
  int j, kmo;
  Edge e;
  Vertex node3, nnodes, taild;
  int nstats;
  double tailattr, headattr;

  nstats  = (int)N_CHANGE_STATS;
  nnodes = N_NODES;
  kmo = (int)INPUT_PARAM[0] - 1;

  /* *** don't forget tail -> head */
    /* edgestate is 1 if edge exists and will disappear
    edgestate is 0 if edge DNE and will appear */
    tailattr = INPUT_ATTRIB[tail-1];
    headattr = INPUT_ATTRIB[head-1];
    taild = -(int)edgestate; /* if edge exists set to -1 because it will be recounted */

    STEP_THROUGH_OUTEDGES(tail, e, node3) { /* step through outedges of tail */
      if(headattr == INPUT_ATTRIB[node3-1]){++taild;}
    }
    for(j=0; j < N_CHANGE_STATS; j++) {
      if (INPUT_ATTRIB[nnodes+j] == tailattr &&
      INPUT_ATTRIB[nnodes+nstats+j] == headattr) {
        change = CHOOSE(taild, kmo);
        CHANGE_STAT[j] += (edgestate ? - change : change);
      }
    }
}

/*****************
 changestat: d_b1starmixhomophily
*****************/
C_CHANGESTAT_FN(c_b1starmixhomophily) {
  double change;
  int j, kmo;
  Edge e;
  Vertex node3, nnodes, taild;
  double tailattr, headattr;

  nnodes = N_NODES;
  kmo = (int)INPUT_PARAM[0] - 1;


  /* *** don't forget tail -> head */
    /* edgestate is 1 if edge exists and will disappear
    edgestate is 0 if edge DNE and will appear */
    tailattr = INPUT_ATTRIB[tail-1];
    headattr = INPUT_ATTRIB[head-1];
    taild = -(int)edgestate; /* if edge exists set to -1 because it will be recounted */

    STEP_THROUGH_OUTEDGES(tail, e, node3) { /* step through outedges of tail */
      if(headattr == INPUT_ATTRIB[node3-1]){++taild;}
    }
    for(j=0; j < N_CHANGE_STATS; j++) {
      if (INPUT_ATTRIB[nnodes+j] == tailattr) {
        change = CHOOSE(taild, kmo);
        CHANGE_STAT[j] += (edgestate ? - change : change);
      }
    }
}

/*****************
 changestat: d_b1twostar
*****************/
C_CHANGESTAT_FN(c_b1twostar) {
  double change;
  int j;
  Edge e;
  Vertex node3, nnodes;
  int nstats;
  double tailattr, headattr, n3attr;

  nstats  = (int)N_CHANGE_STATS;
  nnodes = N_NODES;

  /* *** don't forget tail -> head */
    change = IS_OUTEDGE(tail = tail, head)? -1.0 : 1.0 ;
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
}

/*****************
 changestat: d_b2concurrent
*****************/
C_CHANGESTAT_FN(c_b2concurrent) {
  int echange;
  Vertex b2, b2deg;

  /* *** don't forget tail -> head */
    b2 = head;
    echange = IS_OUTEDGE(tail, b2) ? -1 : 1;
    b2deg = IN_DEG[b2];
    CHANGE_STAT[0] += (b2deg + echange > 1) - (b2deg > 1);
}

/*****************
 changestat: d_b2concurrent_by_attr
*****************/
C_CHANGESTAT_FN(c_b2concurrent_by_attr) {
  /*The inputparams are assumed to be set up as follows:
    The first 2*nstats values are in pairs:  (degree, attrvalue)
    The values following the first 2*nstats values are the nodal attributes.*/
  int j, echange, b2attr;
  Vertex b2, b2deg;


  /* *** don't forget tail -> head */
    b2 = head;
    echange = IS_OUTEDGE(tail, b2) ? -1 : 1;
    b2deg = IN_DEG[b2];
    b2attr = INPUT_PARAM[N_CHANGE_STATS + b2 - 1 - BIPARTITE];
    for(j = 0; j < N_CHANGE_STATS; j++) {
      if (b2attr == INPUT_PARAM[j]) { /* we have attr match */
        CHANGE_STAT[j] += (b2deg + echange > 1) - (b2deg > 1);
      }
    }
}

/*****************
 changestat: d_b2nodematch
*****************/
C_CHANGESTAT_FN(c_b2nodematch) {

  Vertex node3, node4, ninputs;
  int count, exponenttype, matchval, b1attrsize, attrval1, attrval2, diffstatus;
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

    matchval = INPUT_PARAM[head + ninputs - BIPARTITE - 1];
    /* Now count the neighbors of tail whose attribute value equals matchval */
    /* All neighbors of tail are outedges because this is a bipartite network */
    count=0;
    change = 0.0;

    /* RPB */
    /*  double CHANGE[b1attrsize]; */


  if(b1attrsize == 0){

    STEP_THROUGH_OUTEDGES(tail, e, node3) {
      if (INPUT_PARAM[node3 + ninputs - BIPARTITE - 1] == matchval && head != node3) { /* match! */
        ++count;

	// Rprintf("Matching twostar found! %d and %d connect to %d\n==================\n", tail, node3, head);
        if (exponenttype == AlphaType) {
          /* calculate alpha change stat instead of beta change stat. */
          /* Look for number of two-paths connecting head and node3 */
          count = 0;

      STEP_THROUGH_INEDGES(head, e2, node4) {
            // Rprintf("node3=%d, node4=%d, alpha=%f\n", node3,node4,alpha);
            if (node4 != tail) {
              count += IS_OUTEDGE(node4, node3); /* add 1 if node4 connects node3 with head */
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
	    CHANGE_STAT[0] += edgestate ? -change : change;
	    else
	    CHANGE_STAT[matchval-1] += edgestate ? -change : change;

    } else { /* diff=F */
	    CHANGE_STAT[0] += edgestate ? -change : change;
    }
  } else {

 attrval1 = INPUT_PARAM[tail + ninputs + N_NODES + b1attrsize - BIPARTITE - 1];

      STEP_THROUGH_OUTEDGES(tail, e, node3) {

	if (INPUT_PARAM[node3 + ninputs - BIPARTITE - 1] == matchval && head != node3) { /* match! */

	  ++count;

	  // Rprintf("Matching twostar found! %d and %d connect to %d\n==================\n", tail, node3, head);
	  if (exponenttype == AlphaType) {
	    /* calculate alpha change stat instead of beta change stat. */
	    /* Look for number of two-paths connecting tail and node3, not via head */

	    count = 0;

	    STEP_THROUGH_INEDGES(head, e2, node4) {
	      // Rprintf("node3=%d, node4=%d, alpha=%f\n", node3,node4,alpha);
	      if (node4 != tail) {
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
	      CHANGE_STAT[b1attrsize*(matchval-1) + attrval1 - 1] += edgestate ? -change : change;
      } else{
          CHANGE_STAT[attrval1 - 1] += edgestate ? -change : change;
      }

     }
}

/*****************
 changestat: d_b2starmix
*****************/
C_CHANGESTAT_FN(c_b2starmix) {
  double change;
  int j, kmo;
  Edge e;
  Vertex node3, nnodes, headd;
  int nstats;
  double tailattr, headattr;

  nstats  = (int)N_CHANGE_STATS;
  nnodes = N_NODES;
  kmo = (int)INPUT_PARAM[0] - 1;


  /* *** don't forget tail -> head */
    /* edgestate is 1 if edge exists and will disappear
    edgestate is 0 if edge DNE and will appear */
    tailattr = INPUT_ATTRIB[tail-1];
    headattr = INPUT_ATTRIB[head-1];
    headd = -(int)edgestate; /* if edge exists set to -1 because it will be recounted */

    STEP_THROUGH_INEDGES(head, e, node3) { /* step through inedges of head */
      if(tailattr == INPUT_ATTRIB[node3-1]){++headd;}
    }
    for(j=0; j < N_CHANGE_STATS; j++) {
      if (INPUT_ATTRIB[nnodes+j] == tailattr &&
      INPUT_ATTRIB[nnodes+nstats+j] == headattr) {
        change = CHOOSE(headd, kmo);
        CHANGE_STAT[j] += (edgestate ? - change : change);
      }
    }
}

/*****************
 changestat: d_b2starmixhomophily
*****************/
C_CHANGESTAT_FN(c_b2starmixhomophily) {
  double change;
  int j, kmo;
  Edge e;
  Vertex node3, nnodes, headd;
  double tailattr, headattr;

  nnodes = N_NODES;
  kmo = (int)INPUT_PARAM[0] - 1;


  /* *** don't forget tail -> head */
    /* edgestate is 1 if edge exists and will disappear
    edgestate is 0 if edge DNE and will appear */
    tailattr = INPUT_ATTRIB[tail-1];
    headattr = INPUT_ATTRIB[head-1];
    headd = -(int)edgestate; /* if edge exists set to -1 because it will be recounted */

    STEP_THROUGH_INEDGES(head, e, node3) { /* step through inedges of head */
      if(tailattr == INPUT_ATTRIB[node3-1]){++headd;}
    }
    for(j=0; j < N_CHANGE_STATS; j++) {
      if (INPUT_ATTRIB[nnodes+j] == headattr) {
        change = CHOOSE(headd, kmo);
        CHANGE_STAT[j] += (edgestate ? - change : change);
      }
    }
}

/*****************
 changestat: d_b2twostar
*****************/
C_CHANGESTAT_FN(c_b2twostar) {
  double change;
  int j;
  Edge e;
  Vertex node3, nnodes;
  int nstats;
  double tailattr, headattr, n3attr;

  nstats  = (int)N_CHANGE_STATS;
  nnodes = N_NODES;


  /* *** don't forget tail -> head */
    change = edgestate? -1.0 : 1.0 ;
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
}

/*****************
 changestat: d_balance
*****************/
C_CHANGESTAT_FN(c_balance) {
  int a, b, c, d, e, edgecount, t300,
  t210, t120C, t120U, t120D, t201, t030C, t030T, t111U,
  t111D, t021C, t021U, t021D, t102, t012; /* , t003; */
  Vertex node3;


  /* *** don't forget tail -> head */
  if (DIRECTED) { /* directed version */
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
      CHANGE_STAT[0] += edgestate ? -(double)b : (double)b;

  /* *** don't forget tail -> head */
  }else{ /*  undirected */
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
      CHANGE_STAT[0] += edgestate ? -(double)b : (double)b;
  }
}

/*****************
 changestat: d_boundeddegree
*****************/
C_CHANGESTAT_FN(c_boundeddegree) {
  int j, echange;
  Vertex taild, headd=0, deg;
  int nstats = (int)N_CHANGE_STATS;
  Vertex bound = (Vertex)INPUT_PARAM[nstats-1];

  /* *** don't forget tail -> head */
    echange = edgestate ? -1 : 1;
    taild = DEG(tail);
    headd = DEG(head);
    for(j = 0; j+1 < nstats; j++)	{
      deg = (Vertex)INPUT_PARAM[j];
      CHANGE_STAT[j] += (taild + echange == deg) - (taild == deg);
      CHANGE_STAT[j] += (headd + echange == deg) - (headd == deg);
    }
    CHANGE_STAT[nstats-1] += (taild + echange >= bound) - (taild >= bound);
    CHANGE_STAT[nstats-1] += (headd + echange >= bound) - (headd >= bound);
}

/*****************
 changestat: d_boundedidegree
*****************/
C_CHANGESTAT_FN(c_boundedidegree) {
  int j, echange;
  Vertex taild=0, deg;
  int nstats = (int)N_CHANGE_STATS;
  Vertex bound = (Vertex)INPUT_PARAM[nstats-1];

  /* *** don't forget tail -> head */
    echange = edgestate ? -1 : 1;
    taild = IN_DEG[tail];
    for(j = 0; j < N_CHANGE_STATS; j++) {
      deg = (Vertex)INPUT_PARAM[j];
      CHANGE_STAT[j] += (taild + echange == deg) - (taild == deg);
    }
    CHANGE_STAT[nstats-1] += (taild + echange >= bound) - (taild >= bound);
}

/*****************
 changestat: d_boundedistar
*****************/
C_CHANGESTAT_FN(c_boundedistar) {
  double change, headod;
  double newheadod;
  int j, k, bound;
  int p = N_CHANGE_STATS;

  /* *** don't forget tail -> head */
    /* is there an edge for this toggle */
    headod = IN_DEG[head];
    newheadod = headod + (edgestate ? -1 : 1);
    for(j=0; j < p; j++) {
      k =  ((int)INPUT_PARAM[j]);
      bound = (int)INPUT_PARAM[j+p];
      change = MIN(bound,CHOOSE(newheadod, k))-MIN(bound,CHOOSE(headod, k));
      CHANGE_STAT[j] += change;
    }
}

/*****************
 changestat: d_boundedkstar
*****************/
C_CHANGESTAT_FN(c_boundedkstar) {
  double change, tailod, headod;
  double newtailod, newheadod;
  int j, k, bound;
  int p = N_CHANGE_STATS;

  /* *** don't forget tail -> head */
    /* is there an edge for this toggle */
    tailod = DEG(tail);
    newtailod = tailod + (edgestate ? -1 : 1);
    headod = DEG(head);
    newheadod = headod + (edgestate ? -1 : 1);
    for(j=0; j < p; j++) {
      k =  ((int)INPUT_PARAM[j]);
      bound = (int)INPUT_PARAM[j+p];
      change = (MIN(bound,CHOOSE(newtailod, k))-MIN(bound,CHOOSE(tailod, k))) +
      (MIN(bound,CHOOSE(newheadod, k))-MIN(bound,CHOOSE(headod, k)));

      CHANGE_STAT[j] += change; /* (edgestate ? - change : change); */
    }
}

/*****************
 changestat: d_boundedodegree
*****************/
C_CHANGESTAT_FN(c_boundedodegree) {
  int j, echange;
  Vertex taild=0, deg;
  int nstats = (int)N_CHANGE_STATS;
  Vertex bound = (Vertex)INPUT_PARAM[nstats-1];

  /* *** don't forget tail -> head */
    echange = edgestate ? -1 : 1;
    taild = OUT_DEG[tail];
    for(j = 0; j < N_CHANGE_STATS; j++)  {
      deg = (Vertex)INPUT_PARAM[j];
      CHANGE_STAT[j] += (taild + echange == deg) - (taild == deg);
    }
    CHANGE_STAT[nstats-1] += (taild + echange >= bound) - (taild >= bound);
}

/*****************
 changestat: d_boundedostar
*****************/
C_CHANGESTAT_FN(c_boundedostar) {
  double change, tailod;
  double newtailod;
  int j, k, bound;
  int p = N_CHANGE_STATS;

  /* *** don't forget tail -> head */
    /* is there an edge for this toggle */
    tailod = OUT_DEG[tail];
    newtailod = tailod + (edgestate ? -1 : 1);
      for(j=0; j < p; j++) {
        k =  ((int)INPUT_PARAM[j]);
        bound = (int)INPUT_PARAM[j+p];
        change = MIN(bound,CHOOSE(newtailod, k))-MIN(bound,CHOOSE(tailod, k));
        CHANGE_STAT[j] += change;
      }
  }

/*****************
 changestat: d_boundedtriangle
*****************/
Vertex CountTriangles (Vertex tail, Vertex head, int outcount,
                       int incount, Network *nwp);
C_CHANGESTAT_FN(c_boundedtriangle) {
  Edge e;
  Vertex node3;
  double boundedchange, htcount;
  Vertex tailtri, headtri;
  int bound = (int)INPUT_PARAM[0];

  /* *** don't forget tail -> head */
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
    boundedchange = (MIN(headtri+(edgestate ? -1:1)*htcount,bound)-MIN(headtri,bound)+
                    MIN(tailtri+(edgestate ? -1:1)*htcount,bound)-MIN(tailtri,bound));
    CHANGE_STAT[0] += boundedchange;
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
C_CHANGESTAT_FN(c_concurrent) {
  int echange;
  Vertex taildeg, headdeg;

  /* *** don't forget tail -> head */
    echange = edgestate ? -1 : 1;
    taildeg = OUT_DEG[tail];
    headdeg = IN_DEG[head];
    if(!DIRECTED){
      taildeg += IN_DEG[tail];
      headdeg += OUT_DEG[head];
    }
    CHANGE_STAT[0] += (taildeg + echange > 1) - (taildeg > 1);
    CHANGE_STAT[0] += (headdeg + echange > 1) - (headdeg > 1);
}

/*****************
 changestat: d_concurrent_by_attr
*****************/
C_CHANGESTAT_FN(c_concurrent_by_attr) {
  /*The inputparams are assumed to be set up as follows:
    The first 2*nstats values are in pairs:  (degree, attrvalue)
    The values following the first 2*nstats values are the nodal attributes.
  */
  int j, echange, tailattr, headattr;
  Vertex taildeg, headdeg;

  /* *** don't forget tail -> head */
    echange = edgestate ? -1 : 1;
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
}

/*****************
 changestat: d_ctriple
*****************/
C_CHANGESTAT_FN(c_ctriple) {
  Edge e;
  Vertex change, node3;
  int j;
  double tailattr, edgemult;

  /* *** don't forget tail -> head */
    edgemult = edgestate ? -1.0 : 1.0;
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
}

/********************  changestats:  D    ***********/
/*****************
 changestat: d_degcor
*****************/
C_CHANGESTAT_FN(c_degcor) {
  int  echange;
  Vertex taildeg, headdeg, node3;
  Edge e;
  double sigma2;

  sigma2 = INPUT_PARAM[0];
// Rprintf("sigma2 %f\n",sigma2);
    taildeg = DEG(tail);
    headdeg = DEG(head);
    echange = edgestate ? -1 : 1;
    if(echange==1){
     CHANGE_STAT[0] += (taildeg + 1.0)*(headdeg + 1.0);
     STEP_THROUGH_OUTEDGES(head, e, node3) { /* step through outedges of head */
       CHANGE_STAT[0] += DEG(node3);
     }
     STEP_THROUGH_INEDGES(head, e, node3) { /* step through inedges of head */
       CHANGE_STAT[0] += DEG(node3);
     }
     STEP_THROUGH_OUTEDGES(tail, e, node3) { /* step through outedges of tail */
       CHANGE_STAT[0] += DEG(node3);
     }
     STEP_THROUGH_INEDGES(tail, e, node3) { /* step through inedges of tail */
       CHANGE_STAT[0] += DEG(node3);
     }
    }else{
     CHANGE_STAT[0] -= (taildeg)*(headdeg);
     STEP_THROUGH_OUTEDGES(head, e, node3) { /* step through outedges of head */
      if(node3!=tail) CHANGE_STAT[0] -= DEG(node3);
     }
     STEP_THROUGH_INEDGES(head, e, node3) { /* step through inedges of head */
      if(node3!=tail) CHANGE_STAT[0] -= DEG(node3);
     }
     STEP_THROUGH_OUTEDGES(tail, e, node3) { /* step through outedges of tail */
      if(node3!=head) CHANGE_STAT[0] -= DEG(node3);
     }
     STEP_THROUGH_INEDGES(tail, e, node3) { /* step through inedges of tail */
      if(node3!=head) CHANGE_STAT[0] -= DEG(node3);
     }
    }
  CHANGE_STAT[0] *= (2.0/sigma2);
}
S_CHANGESTAT_FN(s_degcor) {
  Vertex taildeg, headdeg;
  Edge e;
  double mu, mu2, sigma2, cross;

  mu = 0.0;
  mu2 = 0.0;
  cross = 0.0;
  for(Vertex tail=1; tail <= N_NODES; tail++) {
    Vertex head;
    STEP_THROUGH_OUTEDGES(tail, e, head) { /* step through outedges of tail */
    taildeg = DEG(tail);
    headdeg = DEG(head);
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
C_CHANGESTAT_FN(c_degcrossprod) {
  int  echange;
  Vertex taildeg, headdeg, node3;
  Edge e;
  double nedges;

  nedges = INPUT_PARAM[0];

    echange = edgestate ? -1 : 1;
    taildeg = DEG(tail);
    headdeg = DEG(head);
    if(echange==1){
     CHANGE_STAT[0] += (taildeg + 1)*(headdeg + 1);
     STEP_THROUGH_OUTEDGES(head, e, node3) { /* step through outedges of head */
      if(node3!=tail) CHANGE_STAT[0] += DEG(node3);
     }
     STEP_THROUGH_INEDGES(head, e, node3) { /* step through inedges of head */
      if(node3!=tail) CHANGE_STAT[0] += DEG(node3);
     }
     STEP_THROUGH_OUTEDGES(tail, e, node3) { /* step through outedges of tail */
      if(node3!=head) CHANGE_STAT[0] += DEG(node3);
     }
     STEP_THROUGH_INEDGES(tail, e, node3) { /* step through inedges of tail */
      if(node3!=head) CHANGE_STAT[0] += DEG(node3);
     }
    }else{
     CHANGE_STAT[0] -= (taildeg)*(headdeg);
     STEP_THROUGH_OUTEDGES(head, e, node3) { /* step through outedges of head */
      if(node3!=tail) CHANGE_STAT[0] -= DEG(node3);
     }
     STEP_THROUGH_INEDGES(head, e, node3) { /* step through inedges of head */
      if(node3!=tail) CHANGE_STAT[0] -= DEG(node3);
     }
     STEP_THROUGH_OUTEDGES(tail, e, node3) { /* step through outedges of tail */
      if(node3!=head) CHANGE_STAT[0] -= DEG(node3);
     }
     STEP_THROUGH_INEDGES(tail, e, node3) { /* step through inedges of tail */
      if(node3!=head) CHANGE_STAT[0] -= DEG(node3);
     }
    }
// Rprintf("N_EDGES %d nedges %f \n",N_EDGES, nedges);
  CHANGE_STAT[0] /= nedges;
}

// A macro indicating whether x is in [from,to)
#define FROM_TO(x, from, to) ((x)>=(from) && (x)<(to))

/*****************
 changestat: d_degrange
*****************/
C_CHANGESTAT_FN(c_degrange) {
  int j, echange;

  /* *** don't forget tail -> head */
      echange=edgestate ? -1:+1;
    Vertex taildeg = DEG(tail), headdeg = DEG(head);
    for(j = 0; j < N_CHANGE_STATS; j++) {
      Vertex from = INPUT_PARAM[2*j], to = INPUT_PARAM[2*j+1];
      CHANGE_STAT[j] += FROM_TO(taildeg + echange, from, to) - FROM_TO(taildeg, from, to);
      CHANGE_STAT[j] += FROM_TO(headdeg + echange, from, to) - FROM_TO(headdeg, from, to);
    }
}

/*****************
 changestat: d_degrange_by_attr
*****************/
C_CHANGESTAT_FN(c_degrange_by_attr) {
  /* The inputparams are assumed to be set up as follows:
  The first 3*nstats values are in triples:  (from, to, attrvalue)
  The values following the first 3*nstats values are the nodal attributes.
  */
  int j;

  /* *** don't forget tail -> head */
      int echange = edgestate ? -1:1;
    Vertex taildeg = DEG(tail), headdeg = DEG(head);
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
}

/*****************
 changestat: d_degrange_w_homophily
*****************/
C_CHANGESTAT_FN(c_degrange_w_homophily) {
  /*  The inputparams are assumed to be set up as follows:
  The first 2*nstats values are the values of degrange
  The values following the first 2*nstats values are the nodal attributes.
  */
  int j;
  Vertex taildeg, headdeg, v;
  double *nodeattr;
  Edge e;

  nodeattr = mtp->inputparams + N_CHANGE_STATS*2 - 1;

  /* *** don't forget tail -> head */
      int tailattr = nodeattr[tail], headattr = nodeattr[head];
    if (tailattr == headattr) { /* They match; otherwise don't bother */
      int echange = edgestate ? -1:1;
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
}

#undef FROM_TO

/*****************
 changestat: d_degree
*****************/
C_CHANGESTAT_FN(c_degree) {
  int j, echange;
  Vertex taildeg, headdeg, deg;

  /* *** don't forget tail -> head */
    echange=edgestate ? -1:+1;
    taildeg = DEG(tail);
    headdeg = DEG(head);
    for(j = 0; j < N_CHANGE_STATS; j++) {
      deg = (Vertex)INPUT_PARAM[j];
      CHANGE_STAT[j] += (taildeg + echange == deg) - (taildeg == deg);
      CHANGE_STAT[j] += (headdeg + echange == deg) - (headdeg == deg);
    }
}

/*****************
 changestat: c_degdist
*****************/
C_CHANGESTAT_FN(c_degdist) {
  int echange = edgestate ? -1:+1;

  Vertex otd = DEG(tail), ohd = DEG(head),
    ntd = otd + echange, nhd = ohd + echange;

  if(ntd > N_CHANGE_STATS || nhd > N_CHANGE_STATS) cutoff_error(mtp);

  if(otd) CHANGE_STAT[otd-1]--;
  if(ohd) CHANGE_STAT[ohd-1]--;
  if(ntd) CHANGE_STAT[ntd-1]++;
  if(nhd) CHANGE_STAT[nhd-1]++;
}

/*****************
 changestat: d_degreepopularity
*****************/
C_CHANGESTAT_FN(c_degreepopularity) {
  double change;

  /* *** don't forget tail -> head */
  change = 0.0;
    Vertex tdeg = DEG(tail);
    Vertex hdeg = DEG(head);
    if(edgestate){
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
  CHANGE_STAT[0]=change;
}

/*****************
 changestat: d_degree_by_attr
*****************/
C_CHANGESTAT_FN(c_degree_by_attr) {
  /* The inputparams are assumed to be set up as follows:
  The first 2*nstats values are in pairs:  (degree, attrvalue)
  The values following the first 2*nstats values are the nodal attributes.
  */
  int j, echange, tailattr, headattr, testattr;
  Vertex taildeg, headdeg, d;

  /* *** don't forget tail -> head */
    echange = edgestate ? -1:1;
    taildeg = DEG(tail);
    headdeg = DEG(head);
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
}

/*****************
 changestat: d_degree_w_homophily
*****************/
C_CHANGESTAT_FN(c_degree_w_homophily) {
  /*  The inputparams are assumed to be set up as follows:
  The first nstats values are the values of degree
  The values following the first nstats values are the nodal attributes.
  */
  int j, echange, tailattr, headattr;
  Vertex taildeg, headdeg, deg, v;
  double *nodeattr;
  Edge e;

  nodeattr = mtp->inputparams + N_CHANGE_STATS - 1;

  /* *** don't forget tail -> head */
    tailattr = (int)nodeattr[tail];
    headattr = (int)nodeattr[head];
    if (tailattr == headattr) { /* They match; otherwise don't bother */
      echange = edgestate ? -1:1;
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
}

/*****************
 changestat: d_dyadcov
*****************/
C_CHANGESTAT_FN(c_dyadcov) {
  double val;
  int refedgestate;
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
      /*Get the initial state of the edge and its reflection*/
      refedgestate = (IS_OUTEDGE(head, tail));

    /* *** don't forget tail -> head */

      /*Get the dyadic covariate*/
      /*    val = INPUT_ATTRIB[(head-1-nrow)+(tail-1)*ncols]; */
      index = (head-1-noffset)*nrow+(tail-1);
      if(index >= 0 && index <= nrow*nrow){
        val = INPUT_ATTRIB[(head-1-noffset)*nrow+(tail-1)];
        /*  Rprintf("tail %d head %d nrow %d ncols %d val %f\n",tail, head, nrow, ncols, val); */

        /*Update the change statistics, as appropriate*/
        if(refedgestate){      /* Reflected edge is present */
          if(edgestate){         /* Toggled edge _was_ present */
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
          if(edgestate){         /* Toggled edge _was_ present */
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
}else{
    /* undirected case (including bipartite) */

    /* *** don't forget tail -> head */
        /*Get the initial edge state*/
      /*Get the covariate value*/
      /*    val = INPUT_ATTRIB[(head-1-nrow)+(tail-1)*ncols]; */
      index = (head-1-noffset)*nrow+(tail-1);
      if(index >= 0 && index <= nrow*((long int)(INPUT_PARAM[0]))){
        val = INPUT_ATTRIB[(head-1-noffset)*nrow+(tail-1)];
        /*Update the change statistic, based on the toggle type*/
        /*  Rprintf("tail %d head %d nrow %d noffset %d val %f\n",tail, head, nrow, noffset, val); */
        /*Update the change statistic, based on the toggle type*/
        CHANGE_STAT[0] += edgestate ? -val : val;
      }
  }
}


/********************  changestats:  G    ***********/

/*****************
 changestat: d_gwdegree
*****************/

#define GWD0(d0) ((decay) ? exp(loneexpd*(d0)) : (d0)==0)

C_CHANGESTAT_FN(c_gwdegree) {
  int  echange=0;
  double decay, loneexpd, change;
  Vertex taild, headd=0;

  decay = INPUT_PARAM[0];
  loneexpd = log1mexp(decay);

  /* *** don't forget tail -> head */
  change = 0.0;
  echange = edgestate ? -1:+1;
  taild = DEG(tail) - edgestate;
  headd = DEG(head) - edgestate;
  change += echange*(GWD0(taild) + GWD0(headd));

  CHANGE_STAT[0] = change;

}

/*****************
 changestat: d_gwdegree_by_attr
*****************/
C_CHANGESTAT_FN(c_gwdegree_by_attr) {
  /*The inputparams are assumed to be set up as follows:
    The first value is the decay parameter (as in Hunter et al, JASA 200?)
    The next sequence of values is the nodal attributes, coded as integers
         from 1 through N_CHANGE_STATS
  */
  int  tailattr, headattr, echange=0;
  double decay, loneexpd;
  Vertex taild, headd=0;

  decay = INPUT_PARAM[0];
  loneexpd = log1mexp(decay);

  /* *** don't forget tail -> head */
    echange = edgestate ? -1:+1;
    taild = DEG(tail) - edgestate;
    tailattr = INPUT_PARAM[tail];
    CHANGE_STAT[tailattr-1] += echange*GWD0(taild);

    headd = DEG(head) - edgestate;
    headattr = INPUT_PARAM[head];
    CHANGE_STAT[headattr-1] += echange*GWD0(headd);

}

/*****************
 changestat: d_gwidegree
*****************/
C_CHANGESTAT_FN(c_gwidegree) {
  double decay, loneexpd, change;
  Vertex headd=0;

  decay = INPUT_PARAM[0];
  loneexpd = log1mexp(decay);
  change = 0.0;

  /* *** don't forget tail -> head */
    headd = IN_DEG[head] - edgestate;
    change += (edgestate? -1.0 : 1.0) * GWD0(headd);
  CHANGE_STAT[0]=change;
}

/*****************
 changestat: d_gwidegree_by_attr
*****************/
C_CHANGESTAT_FN(c_gwidegree_by_attr) {
  /*The inputparams are assumed to be set up as follows:
    The first value is the decay parameter (as in Hunter et al, JASA 2008)
    The next sequence of values is the nodal attributes, coded as integers
         from 1 through N_CHANGE_STATS
  */
  int  headattr, echange;
  double decay, loneexpd;
  Vertex headd;

  decay = INPUT_PARAM[0];
  loneexpd = log1mexp(decay);

  /* *** don't forget tail -> head */
    echange = edgestate ? -1 : 1;
    headd = IN_DEG[head] - edgestate;
    headattr = INPUT_PARAM[head - BIPARTITE]; /* BIPARTITE to make the b2 version a special case. */
    CHANGE_STAT[headattr-1] += echange*GWD0(headd);
}

/*****************
 changestat: d_gwodegree
*****************/
C_CHANGESTAT_FN(c_gwodegree) {
  double decay, loneexpd, change;
  Vertex taild;

  decay = INPUT_PARAM[0];
  loneexpd = log1mexp(decay);
  change = 0.0;

  /* *** don't forget tail -> head */
    taild = OUT_DEG[tail] - edgestate;
    change += (edgestate? -1 : 1) * GWD0(taild);
  CHANGE_STAT[0] = change;
}

/*****************
 changestat: d_gwodegree_by_attr
*****************/
C_CHANGESTAT_FN(c_gwodegree_by_attr) {
  /*The inputparams are assumed to be set up as follows:
    The first value is the decay parameter (as in Hunter et al, JASA 2008)
    The next sequence of values is the nodal attributes, coded as integers
         from 1 through N_CHANGE_STATS
  */
  int  tailattr, echange;
  double decay, loneexpd;
  Vertex taild;

  decay = INPUT_PARAM[0];
  loneexpd = log1mexp(decay);

  /* *** don't forget tail -> head */
    echange = edgestate ? -1 : 1;
    taild = OUT_DEG[tail] - edgestate;
    tailattr = INPUT_PARAM[tail];
    CHANGE_STAT[tailattr-1] += echange*GWD0(taild);
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
C_CHANGESTAT_FN(c_hamming) {
  int  discord;

  Edge wt_net_start= INPUT_PARAM[0]*2+2;
  double defaultval = INPUT_PARAM[wt_net_start-1]; /* Hamming wt for non-edges in cov nwp */
  double *wt_net = INPUT_PARAM+wt_net_start;

  /* *** don't forget tail -> head */

    discord = XOR(dEdgeListSearch(tail, head, INPUT_PARAM), edgestate);

    /* Second, search second network to see if the weight is different from
       defaultval.  In unweighted case, this network is empty. */

    Edge wt_pos = dEdgeListSearch(tail, head, wt_net);
    double val = wt_pos ? wt_net[wt_pos+2*(unsigned int)wt_net[0]] : defaultval;

    CHANGE_STAT[0] += (discord ? -val : val);

}

/********************  changestats:  I    ***********/

// A macro indicating whether x is in [from,to)
#define FROM_TO(x, from, to) ((x)>=(from) && (x)<(to))

/*****************
 changestat: d_idegrange
*****************/
C_CHANGESTAT_FN(c_idegrange) {
  int j, echange;

  /* *** don't forget tail -> head */
      echange=edgestate ? -1:+1;
    Vertex headideg = IN_DEG[head];
    for(j = 0; j < N_CHANGE_STATS; j++) {
      Vertex from = INPUT_PARAM[2*j], to = INPUT_PARAM[2*j+1];
      CHANGE_STAT[j] += FROM_TO(headideg + echange, from, to) - FROM_TO(headideg, from, to);
    }
}

/*****************
 changestat: d_idegrange_by_attr
*****************/
C_CHANGESTAT_FN(c_idegrange_by_attr) {
  /* The inputparams are assumed to be set up as follows:
  The first 3*nstats values are in triples:  (from, to, attrvalue)
  The values following the first 2*nstats values are the nodal attributes.
  */
  int j;

  /* *** don't forget tail -> head */
      int echange = edgestate ? -1:1;
    Vertex headideg = IN_DEG[head];
    int headattr = INPUT_PARAM[3*N_CHANGE_STATS + head - 1 - BIPARTITE];  /* BIPARTITE to make the b2 version a special case. */
    for(j = 0; j < N_CHANGE_STATS; j++){
      Vertex from = INPUT_PARAM[3*j], to = INPUT_PARAM[3*j + 1];
      int testattr = INPUT_PARAM[3*j + 2];
      if (headattr == testattr)  /* we have tail attr match */
        CHANGE_STAT[j] += FROM_TO(headideg + echange, from, to) - FROM_TO(headideg, from, to);
    }
}

/*****************
 changestat: d_idegrange_w_homophily
*****************/
C_CHANGESTAT_FN(c_idegrange_w_homophily) {
  /*  The inputparams are assumed to be set up as follows:
  The first 2*nstats values are the values of idegrange
  The values following the first 2*nstats values are the nodal attributes.
  */
  int j;
  double *nodeattr;
  Edge e;

  nodeattr = mtp->inputparams + N_CHANGE_STATS*2 - 1;

  /* *** don't forget tail -> head */
      int tailattr = nodeattr[tail], headattr = nodeattr[head];
    if (headattr == tailattr) { /* They match; otherwise don't bother */
      int echange = edgestate ? -1:1;
      Vertex headideg=0, v;
      STEP_THROUGH_INEDGES(head, e, v) { headideg += (nodeattr[v]==headattr); }
      for(j = 0; j < N_CHANGE_STATS; j++) {
        Vertex from = INPUT_PARAM[2*j], to = INPUT_PARAM[2*j+1];
        CHANGE_STAT[j] += FROM_TO(headideg + echange, from, to) - FROM_TO(headideg, from, to);
      }
    }
}

#undef FROM_TO

/*****************
 changestat: d_idegree
*****************/
C_CHANGESTAT_FN(c_idegree) {
  int j;

  /* *** don't forget tail -> head */
    int echange = edgestate ? -1 : +1;
    Vertex headd = IN_DEG[head];

    for(j=0; j < N_CHANGE_STATS; j++){
      Vertex deg = INPUT_PARAM[j];
      CHANGE_STAT[j] += (headd + echange == deg) - (headd == deg);
    }
}


/*****************
 changestat: d_idegdist
*****************/
C_CHANGESTAT_FN(c_idegdist) {
  int echange = edgestate ? -1 : +1;
  Vertex ohd = IN_DEG[head],
    nhd = ohd + echange;

  if(nhd > N_CHANGE_STATS) cutoff_error(mtp);

  if(ohd) CHANGE_STAT[ohd-1]--;
  if(nhd) CHANGE_STAT[nhd-1]++;
}


/*****************
 changestat: d_idegree_by_attr
*****************/
C_CHANGESTAT_FN(c_idegree_by_attr) {
  /* The inputparams are assumed to be set up as follows:
  The first 2*nstats values are in pairs:  (degree, attrvalue)
  The values following the first 2*nstats values are the nodal attributes.
  */
  int j, echange, headattr, testattr;
  Vertex headdeg, d;

  /* *** don't forget tail -> head */
    echange=edgestate ? -1 : +1;
    headdeg = IN_DEG[head];
    headattr = INPUT_PARAM[2*N_CHANGE_STATS + head - 1- BIPARTITE];  /* BIPARTITE to make the b2 version a special case. */
    for(j = 0; j < N_CHANGE_STATS; j++) {
      d = (Vertex)INPUT_PARAM[2*j];
      testattr = INPUT_PARAM[2*j + 1];
      if (headattr == testattr)  /* we have head attr match */
        CHANGE_STAT[j] += (headdeg + echange == d) - (headdeg == d);
    }
}

/*****************
 changestat: d_idegree_w_homophily
*****************/
C_CHANGESTAT_FN(c_idegree_w_homophily) {
  /*  The inputparams are assumed to be set up as follows:
  The first nstats values are the values of degree
  The values following the first nstats values are the nodal attributes.
  */
  int j, echange, tailattr, headattr;
  Vertex headdeg, deg, tmp;
  double *nodeattr;
  Edge e;

  nodeattr = mtp->inputparams + N_CHANGE_STATS - 1;

  /* *** don't forget tail -> head */
    tailattr = (int)nodeattr[tail];
    headattr = (int)nodeattr[head];
    if (tailattr == headattr) { /* They match; otherwise don't bother */
      echange=edgestate ? -1 : +1;
      headdeg=0;
      STEP_THROUGH_INEDGES(head, e, tmp){
        headdeg += (nodeattr[tmp]==headattr);
      }
      for(j = 0; j < N_CHANGE_STATS; j++) {
        deg = (Vertex)INPUT_PARAM[j];
        CHANGE_STAT[j] += (headdeg + echange == deg) - (headdeg == deg);
      }
    }
}

/*****************
 changestat: d_idegreepopularity
*****************/
C_CHANGESTAT_FN(c_idegreepopularity) {
  double change;
  Vertex  deg=0;

  /* *** don't forget tail -> head */
  change = 0.0;
    deg = (double)(IN_DEG[head]);
    if(edgestate){
      change -= sqrt(deg);
      change += (deg-1.0)*(sqrt(deg-1.0)-sqrt(deg));
    }else{
      change += sqrt(deg+1.0);
      change += deg*(sqrt(deg+1.0)-sqrt(deg));
  }
  CHANGE_STAT[0]=change;
}

/*****************
 changestat: d_intransitive
*****************/
C_CHANGESTAT_FN(c_intransitive) {
  Edge e;
  Vertex node2;
  double change;

  /* *** don't forget tail -> head */
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
    CHANGE_STAT[0] += edgestate ? -change : change;
/*  Rprintf("tail %d head %d edgestate %d change %f\n",tail,head, change); */
}

/*****************
changestat: d_isolatededges
*****************/
D_CHANGESTAT_FN(d_isolatededges) {
  int i, edgestate;
  Vertex tail, head, neighbor, taild, headd;
  Edge e;

  /* *** don't forget tail -> head */
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    // is there an edge tail -> head?
    edgestate = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i));

    taild = DEG(tail);
    headd = DEG(head);

    if(edgestate) { // we are removing an edge

      // if head and tail both have degree one, then
      // we are removing an isolated edge
      if(taild == 1 && headd == 1)
        CHANGE_STAT[0] -= 1;

      // if tail has degree 2 and has a degree one node other than head as a neighbor,
      // then we are making a non-isolated edge into an isolated edge by removing
      // the edge tail -> head
      if(taild == 2) {
        STEP_THROUGH_OUTEDGES(tail, e, neighbor) {
          if(DEG(neighbor) == 1 && neighbor != head)
            CHANGE_STAT[0] += 1;
        }
        STEP_THROUGH_INEDGES(tail, e, neighbor) {
          if(DEG(neighbor) == 1 && neighbor != head)
            CHANGE_STAT[0] += 1;
        }
      }

      // ditto head
      if(headd == 2) {
        STEP_THROUGH_OUTEDGES(head, e, neighbor) {
          if(DEG(neighbor) == 1 && neighbor != tail)
            CHANGE_STAT[0] += 1;
        }
        STEP_THROUGH_INEDGES(head, e, neighbor) {
          if(DEG(neighbor) == 1 && neighbor != tail)
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
          if(DEG(neighbor) == 1)
            CHANGE_STAT[0] -= 1;
        }
        STEP_THROUGH_INEDGES(tail, e, neighbor) {
          if(DEG(neighbor) == 1)
            CHANGE_STAT[0] -= 1;
        }
      }

      // ditto head
      if(headd == 1) {
        STEP_THROUGH_OUTEDGES(head, e, neighbor) {
          if(DEG(neighbor) == 1)
            CHANGE_STAT[0] -= 1;
        }
        STEP_THROUGH_INEDGES(head, e, neighbor) {
          if(DEG(neighbor) == 1)
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
C_CHANGESTAT_FN(c_isolates) {
  int  echange;
  Vertex taild, headd=0;

  /* *** don't forget tail -> head */
      echange = edgestate ? -1:+1;
      taild = DEG(tail);
      headd = DEG(head);
      CHANGE_STAT[0] += (taild + echange == 0) - (taild == 0);
      CHANGE_STAT[0] += (headd + echange == 0) - (headd == 0);

}

S_CHANGESTAT_FN(s_isolates) {
  /* *** don't forget tail -> head */
  CHANGE_STAT[0] = 0.0;
  for(Vertex tail=1; tail <= N_NODES; tail++){
    if(DEG(tail) == 0)
      CHANGE_STAT[0] ++;
  }
}

/*****************
 changestat: d_istar
*****************/
C_CHANGESTAT_FN(c_istar) {
  double change, headd=0.0;
  int j, kmo;
  Edge e;
  Vertex node3;
  int ninputs, nstats;
  double tailattr;

  ninputs = (int)N_INPUT_PARAMS;
  nstats  = (int)N_CHANGE_STATS;

  /* *** don't forget tail -> head */
  if(ninputs>nstats){
    /* match on attributes */
      /* edgestate is 1 if edge exists and will disappear
      edgestate is 0 if edge DNE and will appear */
      tailattr = INPUT_ATTRIB[tail-1];
      if(tailattr == INPUT_ATTRIB[head-1]){
        headd = -(int)edgestate;
        STEP_THROUGH_INEDGES(head, e, node3) {/* step through inedges of head */
          if(tailattr == INPUT_ATTRIB[node3-1]){++headd;}
        }
        for(j=0; j < N_CHANGE_STATS; j++) {
          kmo = ((int)INPUT_PARAM[j]) - 1;
          change = CHOOSE(headd, kmo);
          CHANGE_STAT[j] += (edgestate ? - change : change);
        }
      }
  }else{
      /* edgestate is 1 if edge exists and will disappear
      edgestate is 0 if edge DNE and will appear */
      headd = IN_DEG[head] - edgestate;
      for(j=0; j < N_CHANGE_STATS; j++) {
        kmo = ((int)INPUT_PARAM[j]) - 1;
        change = CHOOSE(headd, kmo);
        CHANGE_STAT[j] += (edgestate ? - change : change);
      }
  }
}

/********************  changestats:  K    ***********/
/*****************
 changestat: d_kstar
*****************/
C_CHANGESTAT_FN(c_kstar) {
  double change, taild, headd=0.0;
  int j, kmo;
  Edge e;
  Vertex node3;
  int ninputs, nstats;
  double tailattr;

  ninputs = (int)N_INPUT_PARAMS;
  nstats  = (int)N_CHANGE_STATS;

  /* *** don't forget tail -> head */
  if(ninputs>nstats){
    /* match on attributes */
      /* edgestate is 1 if edge exists and will disappear
      edgestate is 0 if edge DNE and will appear */
      tailattr = INPUT_ATTRIB[tail-1];
      if(tailattr == INPUT_ATTRIB[head-1]){
        taild = -(int)edgestate;
        STEP_THROUGH_OUTEDGES(tail, e, node3) {
          if(tailattr == INPUT_ATTRIB[node3-1]){++taild;}
        }
        STEP_THROUGH_INEDGES(tail, e, node3) {
          if(tailattr == INPUT_ATTRIB[node3-1]){++taild;}
        }
        headd = -(int)edgestate;
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
          CHANGE_STAT[j] += (edgestate ? - change : change);
        }
      }
  }else{
    /* *** don't forget tail -> head */
      /* edgestate is 1 if edge exists and will disappear
      edgestate is 0 if edge DNE and will appear */
      taild = DEG(tail) - edgestate;
      headd = DEG(head) - edgestate;
      for(j=0; j < N_CHANGE_STATS; j++)
      {
        kmo = ((int)INPUT_PARAM[j]) - 1;
/*        if (kmo==0) {
          change=1;
        } else { */
          change = CHOOSE(taild, kmo) + CHOOSE(headd, kmo);
/*      } uncomment these few lines to define 1-stars as equivalent to
          edges (currently, each edge is counted as two 1-stars) */
        CHANGE_STAT[j] += (edgestate ? - change : change);
      }

  }
}


/********************  changestats:  L    ***********/
/*****************
 changestat: d_localtriangle
*****************/
C_CHANGESTAT_FN(c_localtriangle) {
  Edge e;
  Vertex node3, nmat;
  double change;

  nmat = (Vertex)(INPUT_PARAM[0]);

  /* *** don't forget tail -> head */
      change = 0.0;

      if(INPUT_PARAM[1+(head-1)+(tail-1)*nmat] == 1.0){
        STEP_THROUGH_OUTEDGES(head, e, node3) { /* step through outedges of head */
	    if(INPUT_PARAM[1+(node3-1)+(tail-1)*nmat] == 1.0 &&
	       INPUT_PARAM[1+(node3-1)+(head-1)*nmat] == 1.0 ){
	      if (DIRECTED){
		if (IS_INEDGE(node3,tail) ) ++change;
		if (IS_OUTEDGE(node3,tail)) ++change;
	      }else{
		if (IS_UNDIRECTED_EDGE(node3,tail)) ++change;
	      }
	    }
	  }

        STEP_THROUGH_INEDGES(head, e, node3) { /* step through inedges of head */
	    if(INPUT_PARAM[1+(node3-1)+(tail-1)*nmat] == 1.0 &&
	       INPUT_PARAM[1+(node3-1)+(head-1)*nmat] == 1.0 ){
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

	CHANGE_STAT[0] += edgestate ? - change : change;

      }
}

/********************  changestats:  M    ***********/
/*****************
 changestat: d_m2star
*****************/
C_CHANGESTAT_FN(c_m2star) {
  int tailid, headod, change;
  int backedgestate;


  /* *** don't forget tail -> head */
      /*  edgestate is 1 if the edge from tail to head  */
      /*   exists and will disappear */
      /*  edgestate is 0 if the edge does not exist */
      backedgestate = (IS_OUTEDGE(head, tail));

      tailid = IN_DEG[tail];
      headod = OUT_DEG[head];
      change = tailid + headod - 2*backedgestate;
      CHANGE_STAT[0] += (edgestate ? -change : change);

}

/*****************
 changestat: d_mutual

 (1,1) -> anything = -1
 anything -> (1,1) = +1
*****************/
C_CHANGESTAT_FN(c_mutual) {
  double matchval, change;
  int j, ninputs, noattr;

  ninputs = N_INPUT_PARAMS - N_NODES;
  noattr = (N_INPUT_PARAMS == 0);

  /* *** don't forget tail -> head */
    if (IS_OUTEDGE(head,tail)) { /* otherwise, no change occurs */
      change = edgestate ? -1.0 : 1.0 ;
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
}

/*****************
 changestat: d_mutual_by_attr
*****************/
C_CHANGESTAT_FN(c_mutual_by_attr) {
  double change;
  int j, ninputs;

  ninputs = N_INPUT_PARAMS - N_NODES;

  /* *** don't forget tail -> head */
    if (IS_OUTEDGE(head,tail)) { /* otherwise, no change occurs */
      change = edgestate ? -1.0 : 1.0 ;
      for (j=0; j<ninputs; j++) {
        if (INPUT_PARAM[tail+ninputs-1] == INPUT_PARAM[j]){CHANGE_STAT[j] += change;}
        if (INPUT_PARAM[head+ninputs-1] == INPUT_PARAM[j]){CHANGE_STAT[j] += change;}
      }
    }
}

/********************  changestats:  N    ***********/
/*****************
 changestat: d_nearsimmelian
*****************/
C_CHANGESTAT_FN(c_nearsimmelian) {
  Vertex node3;
  double change;
  int edgestateth, sc;

  /* *** don't forget tail -> head */
  edgestateth = (!IS_OUTEDGE(head,tail));

  for(node3=1;node3<=N_NODES;node3++){
    if((node3!=tail)&&(node3!=head)){
     sc = edgestateth + (!IS_OUTEDGE(node3,tail));
     if(sc < 2){
      sc += (!IS_OUTEDGE(tail,node3));
      if(sc < 2){
       sc += (!IS_OUTEDGE(node3,head));
       if(sc < 2){
        sc += (!IS_OUTEDGE(head,node3));
        if(sc < 2){
         change=0.0;
         if (sc == 0 && edgestate == 0 ){--change;}
         if (sc == 0 && edgestate == 1 ){++change;}
         if (sc == 1 && edgestate == 0 ){++change;}
         if (sc == 1 && edgestate == 1 ){--change;}
         CHANGE_STAT[0] += change;
	}
       }
      }
     }
    }
   }

}

/********************  changestats:  O    ***********/

// A macro indicating whether x is in [from,to)
#define FROM_TO(x, from, to) ((x)>=(from) && (x)<(to))

/*****************
 changestat: d_odegrange
*****************/
C_CHANGESTAT_FN(c_odegrange) {
  int j, echange;

  /* *** don't forget tail -> head */
      echange=edgestate ? -1:+1;
    Vertex tailodeg = OUT_DEG[tail];
    for(j = 0; j < N_CHANGE_STATS; j++) {
      Vertex from = INPUT_PARAM[2*j], to = INPUT_PARAM[2*j+1];
      CHANGE_STAT[j] += FROM_TO(tailodeg + echange, from, to) - FROM_TO(tailodeg, from, to);
    }
}

/*****************
 changestat: d_odegrange_by_attr
*****************/
C_CHANGESTAT_FN(c_odegrange_by_attr) {
  /* The inputparams are assumed to be set up as follows:
  The first 3*nstats values are in triples:  (from, to, attrvalue)
  The values following the first 2*nstats values are the nodal attributes.
  */
  int j;

  /* *** don't forget tail -> head */
      int echange = edgestate ? -1:1;
    Vertex tailodeg = OUT_DEG[tail];
    int tailattr = INPUT_PARAM[3*N_CHANGE_STATS + tail - 1];
    for(j = 0; j < N_CHANGE_STATS; j++){
      Vertex from = INPUT_PARAM[3*j], to = INPUT_PARAM[3*j + 1];
      int testattr = INPUT_PARAM[3*j + 2];
      if (tailattr == testattr)  /* we have tail attr match */
        CHANGE_STAT[j] += FROM_TO(tailodeg + echange, from, to) - FROM_TO(tailodeg, from, to);
    }
}

/*****************
 changestat: d_odegrange_w_homophily
*****************/
C_CHANGESTAT_FN(c_odegrange_w_homophily) {
  /*  The inputparams are assumed to be set up as follows:
  The first 2*nstats values are the values of odegrange
  The values following the first 2*nstats values are the nodal attributes.
  */
  int j;
  double *nodeattr;
  Edge e;

  nodeattr = mtp->inputparams + N_CHANGE_STATS*2 - 1;

  /* *** don't forget tail -> head */
      int tailattr = nodeattr[tail], headattr = nodeattr[head];
    if (tailattr == headattr) { /* They match; otherwise don't bother */
      int echange = edgestate ? -1:1;
      Vertex tailodeg=0, v;
      STEP_THROUGH_OUTEDGES(tail, e, v) { tailodeg += (nodeattr[v]==tailattr); }
      for(j = 0; j < N_CHANGE_STATS; j++) {
        Vertex from = INPUT_PARAM[2*j], to = INPUT_PARAM[2*j+1];
        CHANGE_STAT[j] += FROM_TO(tailodeg + echange, from, to) - FROM_TO(tailodeg, from, to);
      }
    }
}

#undef FROM_TO

/*****************
 changestat: d_odegree
*****************/
C_CHANGESTAT_FN(c_odegree) {
  int j;

  /* *** don't forget tail -> head */
      int echange = edgestate ? -1 : 1;
    Vertex taild = OUT_DEG[tail];

    for(j=0; j < N_CHANGE_STATS; j++) {
      Vertex deg = INPUT_PARAM[j];
      CHANGE_STAT[j] = (taild + echange == deg) - (taild == deg);
    }
}


/*****************
 changestat: d_odegdist
*****************/
C_CHANGESTAT_FN(c_odegdist) {
  int echange = edgestate ? -1 : +1;
  Vertex otd = OUT_DEG[tail],
    ntd = otd + echange;

  if(ntd > N_CHANGE_STATS) cutoff_error(mtp);

  if(otd) CHANGE_STAT[otd-1]--;
  if(ntd) CHANGE_STAT[ntd-1]++;
}


/*****************
 changestat: d_odegree_by_attr
*****************/
C_CHANGESTAT_FN(c_odegree_by_attr) {
  /* The inputparams are assumed to be set up as follows:
  The first 2*nstats values are in pairs:  (degree, attrvalue)
  The values following the first 2*nstats values are the nodal attributes.
  */
  int j, echange, tailattr, testattr;
  Vertex taildeg, d;

  /* *** don't forget tail -> head */
    echange=edgestate ? -1 : +1;
    taildeg = OUT_DEG[tail];
    tailattr = INPUT_PARAM[2*N_CHANGE_STATS + tail - 1];
    for(j = 0; j < N_CHANGE_STATS; j++) {
      d = (Vertex)INPUT_PARAM[2*j];
      testattr = INPUT_PARAM[2*j + 1];
      if (tailattr == testattr) { /* we have tail attr match */
        CHANGE_STAT[j] += (taildeg + echange == d) - (taildeg == d);
      }
    }
}

/*****************
 changestat: d_odegree_w_homophily
*****************/
C_CHANGESTAT_FN(c_odegree_w_homophily) {
  /*  The inputparams are assumed to be set up as follows:
  The first nstats values are the values of degree
  The values following the first nstats values are the nodal attributes.
  */
  int j;
  double *nodeattr;
  Edge e;

  nodeattr = mtp->inputparams + N_CHANGE_STATS - 1;

  /* *** don't forget tail -> head */
      int tailattr = nodeattr[tail], headattr = nodeattr[head];
    if (tailattr == headattr) { /* They match; otherwise don't bother */
      int echange=edgestate ? -1 : +1;
      Vertex taildeg=0, tmp;
      STEP_THROUGH_OUTEDGES(tail, e, tmp){
        taildeg += (nodeattr[tmp]==tailattr);
      }
      for(j = 0; j < N_CHANGE_STATS; j++) {
        Vertex deg = INPUT_PARAM[j];
        CHANGE_STAT[j] += (taildeg + echange == deg) - (taildeg == deg);
      }
    }
}

/*****************
 changestat: d_opentriad
*****************/
C_CHANGESTAT_FN(c_opentriad) {

  /* *** don't forget tail -> head */
    Vertex node3;
    Edge change = 0, e;
    /* edgestate is 1 if edge exists and will disappear
       edgestate is 0 if edge DNE and will appear */

    // -3 * triangles

    STEP_THROUGH_OUTEDGES(head, e, node3) { /* step through outedges of head */
      change += IS_UNDIRECTED_EDGE(node3,tail);
    }
    STEP_THROUGH_INEDGES(head, e, node3) { /* step through inedges of head */
      change += IS_UNDIRECTED_EDGE(node3,tail);
    }
    CHANGE_STAT[0] += change * (edgestate ? 3.0 : -3.0);


    // +1 * 2-stars

    Vertex taild = DEG(tail) - edgestate;
    Vertex headd = DEG(head) - edgestate;
    change = taild + headd;
    CHANGE_STAT[0] += (edgestate ?  -change : change);

}

/*****************
 changestat: d_ostar
*****************/
C_CHANGESTAT_FN(c_ostar) {
  double change, headd=0.0;
  int j, kmo;
  Edge e;
  Vertex node3;
  int ninputs, nstats;
  double headattr;

  ninputs = (int)N_INPUT_PARAMS;
  nstats  = (int)N_CHANGE_STATS;

  /* *** don't forget tail -> head */
  if(ninputs>nstats){
    /* match on attributes */
      /* edgestate is 1 if edge exists and will disappear
      edgestate is 0 if edge DNE and will appear */
      headattr = INPUT_ATTRIB[head-1];
      if(headattr == INPUT_ATTRIB[tail-1]){
        headd = -(int)edgestate;
        STEP_THROUGH_OUTEDGES(tail, e, node3) { /* step through outedges of head */
          if(headattr == INPUT_ATTRIB[node3-1]){++headd;}
        }
        for(j=0; j < N_CHANGE_STATS; j++) {
          kmo = ((int)INPUT_PARAM[j]) - 1;
          change = CHOOSE(headd, kmo);
          CHANGE_STAT[j] += (edgestate ? - change : change);
        }
      }
    }else{
      /* edgestate is 1 if edge exists and will disappear
      edgestate is 0 if edge DNE and will appear */
      headd = OUT_DEG[tail] - edgestate;
      for(j=0; j < N_CHANGE_STATS; j++) {
        kmo = ((int)INPUT_PARAM[j]) - 1;
        change = CHOOSE(headd, kmo);
        CHANGE_STAT[j] += (edgestate ? - change : change);
      }
  }
}

/*****************
 changestat: d_odegreepopularity
*****************/
C_CHANGESTAT_FN(c_odegreepopularity) {
  double change;
  Vertex  deg=0;

  /* *** don't forget tail -> head */
  change = 0.0;
    deg = (double)(OUT_DEG[tail]);
    if(edgestate){
      change -= sqrt(deg);
      change += (deg-1.0)*(sqrt(deg-1.0)-sqrt(deg));
    }else{
      change += sqrt(deg+1.0);
      change += deg*(sqrt(deg+1.0)-sqrt(deg));
    }
CHANGE_STAT[0]=change;
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
  Vertex taildeg, headdeg;
  Edge e;
  double mu, mu2, mutail, mutail2, sigma2, sigmatail2, cross;

  mu = 0.0;
  mu2 = 0.0;
  mutail = 0.0;
  mutail2 = 0.0;
  cross = 0.0;
  for(Vertex tail=1; tail <= N_NODES; tail++) {
    Vertex head;
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
  Vertex taildeg, headdeg;
  Edge e;
  double mu, mu2, sigma2, cross;
  Vertex tailrank, headrank;
  Vertex *ndeg=R_Calloc(N_NODES+1, Vertex);

  for(Vertex tail=0; tail <= N_NODES; tail++) { ndeg[tail]=0; }
  for(Vertex tail=0; tail < N_NODES; tail++) {
    Vertex head;
   STEP_THROUGH_OUTEDGES(tail, e, head) { /* step through outedges of tail */
    taildeg = DEG(tail);
    headdeg = DEG(head);
    ndeg[taildeg+1]++;
    ndeg[headdeg+1]++;
   }
  }
for(Vertex tail=1; tail <= N_NODES; tail++) {
    ndeg[tail] += ndeg[tail-1];
}
// Rprintf("tail  %d taildeg[tail] %d \n",tail,ndeg[tail]);}

  mu = 0.0;
  mu2 = 0.0;
  cross = 0.0;
  for(Vertex tail=1; tail <= N_NODES; tail++) {
    Vertex head;
   STEP_THROUGH_OUTEDGES(tail, e, head) { /* step through outedges of tail */
    taildeg = DEG(tail);
    headdeg = DEG(head);
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
  R_Free(ndeg);
}

/*****************
 changestat: d_simmelian
*****************/
C_CHANGESTAT_FN(c_simmelian) {
  Edge e;
  Vertex change, node3;

  /* *** don't forget tail -> head */

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

   CHANGE_STAT[0] += edgestate ? -(double)change : (double)change;
   }

}

/*****************
 changestat: d_simmelianties
*****************/
C_CHANGESTAT_FN(c_simmelianties) {
  Edge e, e2;
  Vertex change, node3, node4, first, htflag;

  /* *** don't forget tail -> head */

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
     CHANGE_STAT[0] += edgestate ? -(double)change : (double)change;
   }
}

/********************  changestats:  T    ***********/

/*****************
 changestat: d_threetrail
*****************/
C_CHANGESTAT_FN(c_threetrail) {
  int j, k, change, dchange[4];
  Edge e;
  Vertex node3;
  /* The four values of dchange represent the four different types of
     directed threetrails oriented so that the middle step is always
     "right" (R).  In order:   RRR, RRL, LRR, LRL
                         i.e., >>>  >><  <>>  <><  */


  /* *** don't forget tail -> head */
    /* Step A: Count threetrails in which tail->head is the middle edge */
    dchange[0] = IN_DEG[tail] * OUT_DEG[head]; /* R then R; may count head->tail->head->tail */
    dchange[1] = IN_DEG[tail] * (IN_DEG[head]-edgestate); /* R then L */
    dchange[2] = (OUT_DEG[tail]-edgestate) * OUT_DEG[head]; /* L then R */
    dchange[3] = (OUT_DEG[tail]-edgestate) * (IN_DEG[head]-edgestate); /* L then L */
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
      dchange[0] -= IS_INEDGE(tail, head) * (1 + 2 * edgestate);
      /* head->tail->head->tail is counted in A whenever IS_INEDGE(tail,head) but
         TT->head->tail->head is only counted in B and C when also edgestate */
      for (j = 0; j < N_INPUT_PARAMS; j++) {
        k = (int) INPUT_PARAM[j];
        CHANGE_STAT[j] += (edgestate ? -dchange[k-1] : dchange[k-1]);
      }
    }
    else { /* Undirected case; don't need head->tail->head->tail correction */
      change = dchange[0] + dchange[1] + dchange[2] + dchange[3];
      CHANGE_STAT[0] += (edgestate ? -change : change);
    }
}

/*****************
 changestat: d_transitive
*****************/
C_CHANGESTAT_FN(c_transitive) {
  Edge e;
  Vertex node2;
  double change;

  /* *** don't forget tail -> head */
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
    CHANGE_STAT[0] += edgestate ? -change : change;
//  Rprintf("tail %d head %d edgestate %d change %f C_S[0]=%f\n", tail, head, change,CHANGE_STAT[0]);
}

C_CHANGESTAT_FN(c_transitiveties) {
  int  echange, ochange;
  int L2th, L2tu, L2uh;
  double cumchange;
  double tailattr;


  /* *** don't forget tail -> head */
    cumchange=0.0;
    L2th=0;
    ochange = GETWT(tail, head) ? -1 : 0;
    echange = 2*ochange + 1;
    if(N_INPUT_PARAMS>0){ /* match on attributes */
      tailattr = INPUT_ATTRIB[tail-1];
      if(tailattr == INPUT_ATTRIB[head-1]){
       /* step through outedges of head  */
       EXEC_THROUGH_OUTEDGES(head,  e,  u, {
         if (GETWT(tail, u) && (tailattr == INPUT_ATTRIB[u-1])){
	   L2tu=ochange;
	   /* step through inedges of u */
	   EXEC_THROUGH_INEDGES(u,  f,  v, {
	     if(GETWT(tail, v) && (tailattr == INPUT_ATTRIB[v-1])){
	       L2tu++;
	       if(L2tu>0) {break;}
	     }
	     });
	   cumchange += (L2tu==0);
         }
       });
       /* step through inedges of head */

       EXEC_THROUGH_INEDGES(head,  e,  u, {
         if (GETWT(tail, u) && (tailattr == INPUT_ATTRIB[u-1])){
	   L2th++;
         }
         if (GETWT(u, tail) && (tailattr == INPUT_ATTRIB[u-1])){
	   L2uh=ochange;
	   /* step through outedges of u */
	   EXEC_THROUGH_OUTEDGES(u,  f,  v, {
	     if(GETWT(v, head) && (tailattr == INPUT_ATTRIB[v-1])){
	       L2uh++;
	       if(L2uh>0) {break;}
	     }
	   });
	   cumchange += (L2uh==0) ;
         }
       });
      }
      }else{ /* no attributes */
    /* step through outedges of head  */
    EXEC_THROUGH_OUTEDGES(head,  e,  u, {
      if (GETWT(tail, u)){
	L2tu=ochange;
	/* step through inedges of u */
	EXEC_THROUGH_INEDGES(u,  f,  v, {
	  if(GETWT(tail, v)){
	    L2tu++;
	    if(L2tu>0) {break;}
	  }
	});
	cumchange += (L2tu==0);
      }
    });
    /* step through inedges of head */

    EXEC_THROUGH_INEDGES(head,  e,  u, {
      if (GETWT(tail, u)){
	L2th++;
      }
      if (GETWT(u, tail)){
	L2uh=ochange;
	/* step through outedges of u */
	EXEC_THROUGH_OUTEDGES(u,  f,  v, {
	  if(GETWT(v, head)){
	    L2uh++;
	    if(L2uh>0) {break;}
	  }
	});
	cumchange += (L2uh==0) ;
      }
    });
    }

    cumchange += (L2th>0) ;
//  Rprintf("L2th %d echange %d cumchange %f tail %d head %d\n", L2th, echange, cumchange,tail,head);
    cumchange  = echange*cumchange;
    (CHANGE_STAT[0]) += cumchange;
}

C_CHANGESTAT_FN(c_cyclicalties) {
  int  echange, ochange;
  int L2th, L2tu, L2uh;
  double cumchange;
  double tailattr;


  /* *** don't forget tail -> head */
    cumchange=0.0;
    L2th=0;
    ochange = GETWT(tail, head) ? -1 : 0;
    echange = 2*ochange + 1;
    if(N_INPUT_PARAMS>0){ /* match on attributes */
      tailattr = INPUT_ATTRIB[tail-1];
      if(tailattr == INPUT_ATTRIB[head-1]){
       /* step through outedges of head  */
       EXEC_THROUGH_OUTEDGES(head,  e,  u, {
         if (GETWT(u, tail) && (tailattr == INPUT_ATTRIB[u-1])){
	   L2tu=ochange;
	   /* step through inedges of u */
	   EXEC_THROUGH_INEDGES(u,  f,  v, {
	     if(GETWT(tail, v) && (tailattr == INPUT_ATTRIB[v-1])){
	       L2tu++;
	       if(L2tu>0) {break;}
	     }
	   });
	   cumchange += (L2tu==0);
         }
       });
       /* step through inedges of head */

       EXEC_THROUGH_OUTEDGES(head,  e,  u, {
         if (GETWT(u, tail) && (tailattr == INPUT_ATTRIB[u-1])){
	   L2th++;
         }
         if (GETWT(u, tail) && (tailattr == INPUT_ATTRIB[u-1])){
	   L2uh=ochange;
	   /* step through outedges of u */
	   EXEC_THROUGH_OUTEDGES(u,  f,  v, {
	     if(GETWT(v, head) && (tailattr == INPUT_ATTRIB[v-1])){
	       L2uh++;
	       if(L2uh>0) {break;}
	     }
	   });
	   cumchange += (L2uh==0) ;
         }
       });
      }
      }else{ /* no attributes */
    /* step through outedges of head  */
    EXEC_THROUGH_OUTEDGES(head,  e,  u, {
      if (GETWT(u, tail)){
	L2tu=ochange;
	/* step through inedges of u */
	EXEC_THROUGH_INEDGES(u,  f,  v, {
	  if(GETWT(tail, v)){
	    L2tu++;
	    if(L2tu>0) {break;}
	  }
	});
	cumchange += (L2tu==0);
      }
    });
    /* step through outedges of head */

    EXEC_THROUGH_OUTEDGES(head,  e,  u, {
      if (GETWT(u, tail)){
	L2th++;
      }
      if (GETWT(u, tail)){
	L2uh=ochange;
	/* step through outedges of u */
	EXEC_THROUGH_OUTEDGES(u,  f,  v, {
	  if(GETWT(v, head)){
	    L2uh++;
	    if(L2uh>0) {break;}
	  }
	});
	cumchange += (L2uh==0) ;
      }
    });
    }

    cumchange += (L2th>0) ;
//  Rprintf("L2th %d echange %d cumchange %f tail %d head %d\n", L2th, echange, cumchange,tail,head);
    cumchange  = echange*cumchange;
    (CHANGE_STAT[0]) += cumchange;
}

/*****************
 changestat: d_triadcensus
*****************/
C_CHANGESTAT_FN(c_triadcensus) {
  int j, a, b, c, d, e, edgecount, t300,
    t210, t120C, t120U, t120D, t201, t030C, t030T, t111U,
    t111D, t021C, t021U, t021D, t102, t012, t003;
  Vertex triadtype, node3;

  /* *** don't forget tail -> head */
  if (DIRECTED) {
    /* directed version */
    t300 = 0;
    t210 = 0;
    t120C = 0;  t120U = 0;   t120D = 0;  t201 = 0;
    t030C = 0;  t030T = 0;   t111U = 0;  t111D = 0;
    t021C = 0;  t021U = 0;   t021D = 0;  t102 = 0;
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
	CHANGE_STAT[j] += edgestate ? -(double)t003 : (double)t003;
	break;
      case 2:   CHANGE_STAT[j] += edgestate ? -(double)t012 : (double)t012;
	break;
      case 3:   CHANGE_STAT[j] += edgestate ? -(double)t102 : (double)t102;
	break;
      case 4:   CHANGE_STAT[j] += edgestate ? -(double)t021D : (double)t021D;
	break;
      case 5:   CHANGE_STAT[j] += edgestate ? -(double)t021U : (double)t021U;
	break;
      case 6:   CHANGE_STAT[j] += edgestate ? -(double)t021C : (double)t021C;
	break;
      case 7:   CHANGE_STAT[j] += edgestate ? -(double)t111D : (double)t111D;
	break;
      case 8:   CHANGE_STAT[j] += edgestate ? -(double)t111U : (double)t111U;
	break;
      case 9:   CHANGE_STAT[j] += edgestate ? -(double)t030T : (double)t030T;
	break;
      case 10:  CHANGE_STAT[j] += edgestate ? -(double)t030C : (double)t030C;
	break;
      case 11:  CHANGE_STAT[j] += edgestate ? -(double)t201 : (double)t201;
	break;
      case 12:  CHANGE_STAT[j] += edgestate ? -(double)t120D : (double)t120D;
	break;
      case 13:  CHANGE_STAT[j] += edgestate ? -(double)t120U : (double)t120U;
	break;
      case 14:  CHANGE_STAT[j] += edgestate ? -(double)t120C : (double)t120C;
	break;
      case 15:  CHANGE_STAT[j] += edgestate ? -(double)t210 : (double)t210;
	break;
      case 16:  CHANGE_STAT[j] += edgestate ? -(double)t300 : (double)t300;
	break;
      }
    }
  } else {
    /*  undirected */

    /* *** don't forget tail -> head */
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
	CHANGE_STAT[j] += edgestate ? -(double)t003 : (double)t003;
	break;
      case 2:  CHANGE_STAT[j] += edgestate ? -(double)t102 : (double)t102;
	break;
      case 3:  CHANGE_STAT[j] += edgestate ? -(double)t201 : (double)t201;
	break;
      case 4:  CHANGE_STAT[j] += edgestate ? -(double)t300 : (double)t300;
	break;
      }
    }
  }
}


/*****************
 changestat: d_tripercent
*****************/
C_CHANGESTAT_FN(c_tripercent) {
  Edge e, e2;
  Vertex node1, node2, node3;
  int j;
  Edge triwith, triwithout;
  Edge degreewith, degreewithout, twostarwith, twostarwithout;
  int ninputs = N_INPUT_PARAMS - N_NODES;
  int MatchingOnAttribute = (ninputs>0);
  double *attr=INPUT_PARAM, ratiowith, ratiowithout;

  if (MatchingOnAttribute)
    attr = INPUT_PARAM + (ninputs-1); /* ptr to vertex attributes */

  /* *** don't forget tail -> head */
  if (!edgestate) TOGGLE(tail, head); /* turn on the edge if it's missing */
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
      CHANGE_STAT[j] += (ratiowith-ratiowithout)*(edgestate? -100.0 : 100.0);
    }
    if (!edgestate) TOGGLE(tail, head); /* reset dyad to original state */
}

/*****************
 changestat: d_ttriple
*****************/
C_CHANGESTAT_FN(c_ttriple) {
  Edge e;
  Vertex change, node3;
  int j;
  double tailattr, edgemult;

  /* *** don't forget tail -> head */
    edgemult = edgestate ? -1.0 : 1.0;
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
}
