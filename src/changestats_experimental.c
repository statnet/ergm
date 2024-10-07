/*  File src/changestats_experimental.c in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2024 Statnet Commons
 */
#include "changestats_experimental.h"


/********************  changestats:  A    ***********/
/*****************
 changestat: d_b1kappa
*****************/
D_CHANGESTAT_FN(d_b1kappa){
  int i, j, echange=0;
  double nedges, change, iar0, far0;
  Vertex tail, head, taild, iak2, fak2, *od;
  Vertex nb1;

  od=OUT_DEG;
  nb1 = BIPARTITE;
  change = 0.0;
  FOR_EACH_TOGGLE(i) {
    echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 1;
    iak2=0;
    for (j=1; j<=nb1; j++) {      
      fak2 = od[j];
      iak2 += fak2*(fak2-1);
    }
    taild = od[tail] + (echange-1)/2;
    fak2 = iak2 + echange*2*taild;
    nedges = (double)(N_EDGES);
    iar0 = (N_EDGES==0) ? 0.0 : (iak2*1.0/nedges);
    far0 = (((N_EDGES)+echange)==0) ? 0.0 : (fak2*1.0/(nedges+echange));
    change += far0 - iar0;
    // Rprintf("tail %d head %d nnodes %d nedges %f iak2 %d fak2 %d iar0 %f far0 %f change %f\n",tail,head, nnodes,  nedges, iak2, fak2, iar0, far0, change);
    // Rprintf("tail %d head %d nnodes %d nedges %f iek2 %d fek2 %d ier0 %f fer0 %f change %f\n",tail,head, nnodes,  nedges, iek2, fek2, ier0, fer0, change);
    TOGGLE_IF_MORE_TO_COME(i);
  }
  CHANGE_STAT[0] = change;
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_altistar
*****************/
D_CHANGESTAT_FN(d_altistar) {
  int i, echange=0;
  double lambda, oneexpl, change;
  Vertex tail, head, headd=0, *id;
  
  id=IN_DEG;
  change = 0.0;
  lambda = INPUT_PARAM[0];
  oneexpl = 1.0-1.0/lambda;
  
  FOR_EACH_TOGGLE(i) {
    echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 1;
    headd = id[head] + (echange - 1)/2;
    if(headd!=0){
      change += echange*(1.0-pow(oneexpl,(double)headd));
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  CHANGE_STAT[0] = change*lambda;
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_altostar
*****************/
D_CHANGESTAT_FN(d_altostar) {
  int i, echange=0;
  double lambda, oneexpl, change;
  Vertex tail, head, taild=0, *od;
  
  od=OUT_DEG;
  change = 0.0;
  lambda = INPUT_PARAM[0];
  oneexpl = 1.0-1.0/lambda;
  
  FOR_EACH_TOGGLE(i) {
    echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 1;
    taild = od[tail] + (echange - 1)/2;
    if(taild!=0){
      change += echange*(1.0-pow(oneexpl,(double)taild));
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  CHANGE_STAT[0] = change*lambda;
  UNDO_PREVIOUS_TOGGLES(i);
}

/********************  changestats:  B    ***********/
/*****************
 changestat: d_berninhom
 Changescores for inhomogeneous Bernoulli graphs
*****************/
D_CHANGESTAT_FN(d_berninhom) {
  int edgestate, i;
  Vertex tail, head, n;

  n=N_NODES;
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    edgestate = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i));
    /*Lower trianglization of adjacency matrix implies ith element corresponds
    to (row-1) + (col-1)*(n-1) - choose(col,2).*/
    /*Rprintf("head=%ld, tail=%ld, cell=%ld, nstats=%d, state=%d\n",head,tail, (head-1)+(tail-1)*(n-1)-tail*(tail-1)/2-1, nstats, edgestate);
    Rprintf("\tdstats content=%f\n",CHANGE_STAT[(head-1)+(tail-1)*(n-1)-tail*(tail-1)/2-1]);*/
    CHANGE_STAT[(head-1)+(tail-1)*(n-1)-tail*(tail-1)/2-1] += edgestate ? - 1 : 1;
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_biduration
*****************/
D_CHANGESTAT_FN(d_biduration)
{
  Vertex tail, head, tailtail, tailhead;
  int i, k, nprevedge, edgestate, discord, lookmore;
  Vertex /* nb2, */ nb1;
  double change=0.0;

  nprevedge = (int)((INPUT_PARAM[0]));
  nb1 = (int)((INPUT_PARAM[1]));
  // nb1 = BIPARTITE;
  // nb2 = N_NODES - BIPARTITE;
  /* nb2 = N_NODES - nb1; */
  
  // Rprintf("nprevedge %d\n", nprevedge);
  // Rprintf("BIPARTITE %d nb2 %f\n", BIPARTITE, BIPARTITE);
  // Rprintf("nb1 %d nb2 %d\n", nb1, nb2);
  // for (k=1; k<=nprevedge; k++) {
    // Rprintf("k %d x0.tail %d x0.head %d\n", k,
    //	      (Vertex)(mtp->attrib[          k]),
    //	      (Vertex)(mtp->attrib[nprevedge+k]));
  // }
  
  CHANGE_STAT[0] = 0.0;
  
  FOR_EACH_TOGGLE(i) {
    /*Get the initial state of the edge and its reflection*/
    edgestate=IS_OUTEDGE(tail=TAIL(i), head=HEAD(i));
    if(tail > head){
      tailtail = head;
      tailhead = tail;
    }else{
      tailtail = tail;
      tailhead = head;
    }
    // 1 = edge removed -1 = edge added
    discord = edgestate ? 1 : -1;
    
    // Rprintf("nprevedge %d\n", nprevedge);    
    k=1;
    lookmore=1;
    while(lookmore && k <= nprevedge){
      // Rprintf("tail %d head %d tailtail %d tailhead %d\n",
      //	      tail,head,tailtail,tailhead);
      // Rprintf("tail %d head %d tailtail %d tailhead %d\n",
      //	      tailhead,tailtail,(Vertex)(mtp->attrib[            k]),
      //	            (Vertex)(mtp->attrib[nprevedge + k]));
      if(tailhead == (Vertex)(mtp->attrib[            k]) &&
        tailtail == (Vertex)(mtp->attrib[nprevedge + k])
      ){
        // Rprintf("tail %d head %d tailtail %d tailhead %d\n",
        //	      tail,head,(Vertex)(mtp->attrib[            k]),
        //	          (Vertex)(mtp->attrib[nprevedge + k]));
        
        /*If the proposed edge existed in x0, get dissolution rate */
        /* lookmore is location of (tailtail, tailhead) in the dissolve matrix */
        // lookmore=2*nprevedge+(tailtail-1)+(tailhead-1-nb2)*nb1+1;
        lookmore=2*nprevedge+(tailtail-1)+(tailhead-1-nb1)*nb1+1;
        /* change is NEGATIVE the dissolve rate if an edge exists */
        // Rprintf("tailtail %d tailhead %d lookmore %d att %f\n",tailtail, tailhead, lookmore, (double)(mtp->attrib[lookmore]));
        // change=-discord*(double)(mtp->attrib[lookmore]);
        change=-discord;
        // Rprintf("lookmore %d change %f discord %d \n",lookmore, change, discord);
        lookmore=0;
        /*Update the change statistics, as appropriate*/
        /*change is the number of edges dissolved*/
        CHANGE_STAT[0] += (double)change;
        // if(change!=0){
          // Rprintf("tail %d head %d tailtail %d tailhead %d edgestate %d beta %f\n",
          //	      tail,head,tailtail,tailhead, edgestate, (double)(mtp->attrib[lookmore]));
          // Rprintf("edgestate %d change %f discord %d\n",
          //	      edgestate, change, discord);
        // }
      }else{
        ++k;
      }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }                              
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_bimix
*****************/
D_CHANGESTAT_FN(d_bimix){

  int matchvaltail, matchvalhead;
  Vertex tail, head, /* nnodes, */ nstats;
  int i, j, edgestate=0;

  nstats = N_CHANGE_STATS;
  /* nnodes = (N_INPUT_PARAMS)-2*nstats; */
  // Rprintf("ninput %d nstats %d\n", N_INPUT_PARAMS, mtp->nstats);
  // Rprintf("nodes %d\n", nnodes);
  // for (tail=0; tail<nnodes; tail++){
    // Rprintf("match tail %d att tail %f\n", 
    //	    tail, mtp->attrib[tail+nstats]);
  // }
  // for (tail=0; tail<nstats; tail++){
    // Rprintf("match tail %d in tail %f in head %f\n", 
    //	    tail, INPUT_PARAM[tail], INPUT_PARAM[tail+nstats]);
  // }
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    tail=TAIL(i);
    head=HEAD(i);
    matchvaltail = mtp->attrib[tail-1+nstats];
    matchvalhead = mtp->attrib[head-1+nstats];
    // if(matchvalhead < matchvaltail){
      // matchswap = matchvaltail;
      // matchvaltail = matchvalhead;
      // matchvalhead = matchswap;
    // }       
    edgestate=IS_OUTEDGE(tail, head);
    for (j=0; j<nstats; j++) {
      if(matchvaltail==INPUT_PARAM[nstats+j] &&
	      matchvalhead==INPUT_PARAM[       j]
      ){CHANGE_STAT[j] += edgestate ? -1.0 : 1.0;}
	  }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_bkappa
*****************/
D_CHANGESTAT_FN(d_bkappa)  {
  int i, j, echange=0;
  double nedges, change, iar0, far0, ier0, fer0;
  Vertex tail, head, taild, headd=0, iak2, fak2, iek2, fek2, nnodes, *id, *od;
  Vertex /* nb2, */ nb1;
  
  id=IN_DEG;
  od=OUT_DEG;
  nnodes = N_NODES;
  nb1 = BIPARTITE;
  /* nb2 = nnodes - nb1; */
  
  change = 0.0;
  FOR_EACH_TOGGLE(i) {
    echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 1;
    iak2=0;
    for (j=1; j<=nb1; j++) {      
      fak2 = od[j];
      iak2 += fak2*(fak2-1);
    }
    iek2=0;
    for (j=nb1+1; j<=nnodes; j++) {      
      fek2 = id[j];
      iek2 += fek2*(fek2-1);
    }
    taild = od[tail] + (echange-1)/2;
    headd = id[head] + (echange-1)/2;
    fak2 = iak2 + echange*2*taild;
    fek2 = iek2 + echange*2*headd;
    nedges = (double)(N_EDGES);
    iar0 = (N_EDGES==0) ? 0.0 : (iak2*1.0/nedges);
    far0 = (((N_EDGES)+echange)==0) ? 0.0 : (fak2*1.0/(nedges+echange));
    ier0 = (N_EDGES==0) ? 0.0 : (iek2*1.0/nedges);
    fer0 = (((N_EDGES)+echange)==0) ? 0.0 : (fek2*1.0/(nedges+echange));
    change += sqrt(far0*fer0) - sqrt(iar0*ier0);
    // Rprintf("tail %d head %d nnodes %d nedges %f iar0 %f far0 %f change %f\n",tail,head, nnodes,  nedges, iar0, far0, change);
    // Rprintf("tail %d head %d nnodes %d nedges %f ier0 %f fer0 %f change %f\n",tail,head, nnodes,  nedges, ier0, fer0, change);
    TOGGLE_IF_MORE_TO_COME(i);
  }
  CHANGE_STAT[0] = change;
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_degreep
*****************/
D_CHANGESTAT_FN(d_degreep) 
{
  int i, j, echange;
  Vertex tail, head, taildeg, headdeg, deg, *id, *od;

  id=IN_DEG;
  od=OUT_DEG;
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    echange=IS_OUTEDGE(tail=TAIL(i), head=HEAD(i))? -1:1;
    taildeg = od[tail] + id[tail];
    headdeg = od[head] + id[head];
    for(j = 0; j < mtp->nstats; j++) {
      deg = (Vertex)INPUT_PARAM[j];
      CHANGE_STAT[j] += ((taildeg + echange == deg) - (taildeg == deg))/(double)N_NODES;
      CHANGE_STAT[j] += ((headdeg + echange == deg) - (headdeg == deg))/(double)N_NODES;
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_degreep_by_attr
*****************/
D_CHANGESTAT_FN(d_degreep_by_attr) 
{
  /* The inputparams are assumed to be set up as follows:
  The first 2*nstats values are in pairs:  (degree, attrvalue)
  The values following the first 2*nstats values are the nodal attributes.
  */
  int i, j, echange, tailattr, headattr, testattr;
  Vertex tail, head, taildeg, headdeg, d, *id, *od;
  
  id=IN_DEG;
  od=OUT_DEG;
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    echange=IS_OUTEDGE(tail=TAIL(i), head=HEAD(i))? -1:1;
    taildeg = od[tail] + id[tail];
    headdeg = od[head] + id[head];
    tailattr = INPUT_PARAM[2*mtp->nstats + tail - 1]; 
    headattr = INPUT_PARAM[2*mtp->nstats + head - 1]; 
    for(j = 0; j < mtp->nstats; j++) {
      d = (Vertex)INPUT_PARAM[2*j];
      testattr = INPUT_PARAM[2*j + 1]; 
      if (tailattr == testattr)  /* we have tail attr match */
        CHANGE_STAT[j] += ((taildeg + echange == d) - (taildeg == d))/(double)N_NODES;
      if (headattr == testattr)  /* we have head attr match */
        CHANGE_STAT[j] += ((headdeg + echange == d) - (headdeg == d))/(double)N_NODES;
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_degreep_w_homophily
*****************/
D_CHANGESTAT_FN(d_degreep_w_homophily) 
{
  /*  The inputparams are assumed to be set up as follows:
  The first nstats values are the values of degree
  The values following the first nstats values are the nodal attributes.
  */
  int i, j, echange, tailattr, headattr;
  Vertex tail, head, taildeg, headdeg, deg, tmp;
  double *nodeattr;
  Edge e;

  nodeattr = INPUT_PARAM + mtp->nstats - 1;  
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    tail=TAIL(i);
    head=HEAD(i);
    tailattr = (int)nodeattr[tail];
    headattr = (int)nodeattr[head];    
    if (tailattr == headattr) { /* They match; otherwise don't bother */
      echange=IS_OUTEDGE(tail, head)? -1:1;
      taildeg=headdeg=0;
      for(e = EdgetreeMinimum(nwp->outedges, tail);
      (tmp = nwp->outedges[e].value) != 0;
      e = EdgetreeSuccessor(nwp->outedges, e)) {
        taildeg += (nodeattr[tmp]==tailattr);
      }
      for(e = EdgetreeMinimum(nwp->inedges, tail);
      (tmp = nwp->inedges[e].value) != 0;
      e = EdgetreeSuccessor(nwp->inedges, e)) {
        taildeg += (nodeattr[tmp]==tailattr);
      }
      for(e = EdgetreeMinimum(nwp->outedges, head);
      (tmp = nwp->outedges[e].value) != 0;
      e = EdgetreeSuccessor(nwp->outedges, e)) {
        headdeg += (nodeattr[tmp]==headattr);
      }
      for(e = EdgetreeMinimum(nwp->inedges, head);
      (tmp = nwp->inedges[e].value) != 0;
      e = EdgetreeSuccessor(nwp->inedges, e)) {
        headdeg += (nodeattr[tmp]==headattr);
      }
      for(j = 0; j < mtp->nstats; j++) {
        deg = (Vertex)INPUT_PARAM[j];
        CHANGE_STAT[j] += ((taildeg + echange == deg) - (taildeg == deg))/(double)N_NODES;
	CHANGE_STAT[j] += ((headdeg + echange == deg) - (headdeg == deg))/(double)N_NODES;
      }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_dissolve
*****************/
D_CHANGESTAT_FN(d_dissolve) {
  int edgestate, i;
  Vertex tail, head;
  
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i) {
      edgestate = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i));
      CHANGE_STAT[0] += edgestate ? - 1 : 1;
    TOGGLE_IF_MORE_TO_COME(i);
    }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_duration
*****************/
D_CHANGESTAT_FN(d_duration) {
  Vertex tail, head, tailtail, tailhead;
  int i, k, ntailedge, edgestate, discord, lookmore;
  int ndyads, nnodes;
  double change=0.0;
  
  ntailedge = (int)((INPUT_PARAM[0]));
  nnodes = N_NODES;
  ndyads = N_NODES * N_NODES;
  
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i) {
    /*Get the initial state of the edge and its reflection*/
    edgestate=IS_OUTEDGE(tail=TAIL(i), head=HEAD(i));
    if(!nwp->directed_flag && tail < head){
      tailtail = head;
      tailhead = tail;
    }else{
      tailtail = tail;
      tailhead = head;
    }
    discord = edgestate ? 1 : -1;
    
    k=0;
    lookmore=1;
    while(lookmore && k < ntailedge){
      if(tailtail == (Vertex)(mtp->attrib[         k]) &&
        tailhead == (Vertex)(mtp->attrib[ntailedge + k])
      ){
        /*Get dissolution rate */
        lookmore=2*ntailedge+ndyads+(tailtail-1)+(tailhead-1)*nnodes;
        change=discord*(double)(mtp->attrib[lookmore]);
        lookmore=0;
      }else{
        ++k;
      }
    }
    if(lookmore){
      /* Get formation rate */
      lookmore=2*ntailedge+(tailtail-1)+(tailhead-1)*nnodes;
      change=-discord*(double)(mtp->attrib[lookmore]);
    }
    /*Update the change statistics, as appropriate*/
    CHANGE_STAT[0] += change;
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/********************  changestats:  C    ***********/

/********************  changestats:  D    ***********/

/********************  changestats:  E    ***********/
/*****************
 changestat: d_b2kappa
*****************/
D_CHANGESTAT_FN(d_b2kappa)  {
  int i, j, echange=0;
  double nedges, change, ier0, fer0;
  Vertex tail, head, headd=0, iek2, fek2, nnodes, *id;
  Vertex nb1;

  id=IN_DEG;
  nnodes = N_NODES;
  nb1 = BIPARTITE;
  
  change = 0.0;
  FOR_EACH_TOGGLE(i) {
    echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : +1;
    iek2=0;
    for (j=nb1+1; j<=nnodes; j++) {      
      fek2 = id[j];
      iek2 += fek2*(fek2-1);
    }
    headd = id[head] + (echange-1)/2;
    fek2 = iek2 + echange*2*headd;
    nedges = (double)(N_EDGES);
    ier0 = (N_EDGES==0) ? 0.0 : (iek2*1.0/nedges);
    fer0 = (((N_EDGES)+echange)==0) ? 0.0 : (fek2*1.0/(nedges+echange));
    change += fer0 - ier0;
    // Rprintf("tail %d t %d nnodes %d nedges %f iak2 %d fak2 %d iar0 %f far0 %f change %f\n",tail,t, nnodes,  nedges, iak2, fak2, iar0, far0, change);
    // Rprintf("tail %d t %d nnodes %d nedges %f iek2 %d fek2 %d ier0 %f fer0 %f change %f\n",tail,t, nnodes,  nedges, iek2, fek2, ier0, fer0, change);
    TOGGLE_IF_MORE_TO_COME(i);
  }
  CHANGE_STAT[0] = change;
  UNDO_PREVIOUS_TOGGLES(i);
}

/********************  changestats:  F    ***********/
/*****************
 changestat: d_factor
*****************/
D_CHANGESTAT_FN(d_factor)  {
  double val, s;
  Vertex tail, head;
  long int nlevels;
  int i, j;
  
  nlevels = (long int)(INPUT_PARAM[0]);

  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    s = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1.0 : 1.0;
    for (j=0; j<(mtp->nstats); j++) {
      /*Get the covariate value*/
      val = INPUT_PARAM[1+j+(tail-1)*nlevels];
      // + INPUT_PARAM[1+j+(head-1)*nlevels];
      CHANGE_STAT[j] += s*val;
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/********************  changestats:  G    ***********/
/*****************
 changestat: d_geodegree
*****************/
D_CHANGESTAT_FN(d_geodegree) {
  int i, echange;
  double alpha;
  Vertex tail, head, taild, headd=0, *id, *od;
                                               
  id=IN_DEG;
  od=OUT_DEG;
  CHANGE_STAT[0] = 0.0;
  alpha = INPUT_PARAM[0];
  
  FOR_EACH_TOGGLE(i) {
    echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 1;
    taild = od[tail] + id[tail] + (echange - 1)/2;
    headd = od[head] + id[head] + (echange - 1)/2;
    CHANGE_STAT[0] += echange*(exp(-alpha*taild)+exp(-alpha*headd));
    TOGGLE_IF_MORE_TO_COME(i);
  }
  CHANGE_STAT[0] = (CHANGE_STAT[0])*(exp(-alpha)-1.0);
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_geospartner
*****************/
D_CHANGESTAT_FN(d_geospartner) {
  Edge e, f;
  int i, echange;
  int L2th, L2tu, L2uh;
  Vertex tail, head, u, v;
  double alpha, cumchange;
  
  CHANGE_STAT[0] = 0.0;
  alpha = INPUT_PARAM[0];
  FOR_EACH_TOGGLE(i) {
    cumchange=0.0;
    L2th=0;
    echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 1;
    STEP_THROUGH_OUTEDGES(head, e, u) {
      if (IS_OUTEDGE(MIN(u,tail), MAX(u,tail))) {
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
        cumchange += exp(-alpha*L2tu)+exp(-alpha*L2uh);
      }
    }
    STEP_THROUGH_INEDGES(head, e, u) {
      if (IS_OUTEDGE(MIN(u,tail), MAX(u,tail))) {
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
        cumchange += exp(-alpha*L2tu)+exp(-alpha*L2uh);
      }
    }
    cumchange  = cumchange*(exp(-alpha*echange)-1.0);
    cumchange += echange*(exp(-alpha*L2th)-1.0);
    CHANGE_STAT[0] -= cumchange;
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_gwb1
*****************/
D_CHANGESTAT_FN(d_gwb1) {
  int i, echange=0;
  double alpha, oneexpa, change;
  Vertex tail, head, taild=0, *od;
  /* int nb1 , nb2 

  nb1 = (int)INPUT_PARAM[0];
  nb2 = (N_NODES) - nb1; */
  
  //id=IN_DEG;
  od=OUT_DEG;
  change = 0.0;
  alpha = INPUT_PARAM[1];
  oneexpa = 1.0-exp(-alpha);
  
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i) {
    echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 1;
    taild = od[tail] + (echange - 1)/2;
    if(taild!=0){
      change += echange*(1.0-pow(oneexpa,(double)taild));
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  CHANGE_STAT[0] = change*exp(alpha);
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_gwd
*****************/
D_CHANGESTAT_FN(d_gwd) {
  int i, echange;
  double alpha;
  Vertex tail, head, taild, headd=0, *id, *od;
  
  id=IN_DEG;
  od=OUT_DEG;
  CHANGE_STAT[0] = 0.0;
  alpha = INPUT_PARAM[0];
  
  FOR_EACH_TOGGLE(i) {
    echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 1;
    taild = od[tail] + id[tail] + (echange - 1)/2;
    headd = od[head] + id[head] + (echange - 1)/2;
    CHANGE_STAT[0] += echange*(exp(-alpha*taild)+exp(-alpha*headd)); 
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_gwdegree706
*****************/
D_CHANGESTAT_FN(d_gwdegree706)  {
  /* Slight modification to the parameterization in d_gwdegree */
  int i, echange=0;
  double decay, oneexpd, change;
  Vertex tail, head, taild, headd=0, *id, *od;
  
  id=IN_DEG;
  od=OUT_DEG;
  decay = INPUT_PARAM[0];
  oneexpd = 1.0-exp(-decay);
  
  change = 0.0;
  FOR_EACH_TOGGLE(i) {
    echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 1;
    change += 4.0*echange;
    taild = od[tail] + id[tail] + (echange - 1)/2;
    headd = od[head] + id[head] + (echange - 1)/2;
    if(taild!=0){
      change -= echange*(1.0-pow(oneexpd,(double)taild));
    }
    if(headd!=0){
      change -= echange*(1.0-pow(oneexpd,(double)headd));
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  CHANGE_STAT[0] = change;
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_gwdegreealpha
*****************/
D_CHANGESTAT_FN(d_gwdegreealpha)  {
  int i, echange=0;
  double alpha, oneexpa, change;
  Vertex tail, head, taild, headd=0, *id, *od;
  
  id=IN_DEG;
  od=OUT_DEG;
  change = 0.0;
  alpha = INPUT_PARAM[0];
  oneexpa = 1.0-exp(-alpha);
  
  FOR_EACH_TOGGLE(i) {
    echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 1;
    taild = od[tail] + id[tail] + (echange - 1)/2;
    headd = od[head] + id[head] + (echange - 1)/2;
    if(taild!=0){
      change += echange*(1.0-pow(oneexpa,(double)taild));
    }
    if(headd!=0){
      change += echange*(1.0-pow(oneexpa,(double)headd));
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  CHANGE_STAT[0] = change*exp(alpha);
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_gwdegreelambda
*****************/
D_CHANGESTAT_FN(d_gwdegreelambda)  {
  int i, echange=0;
  double lambda, oneexpl, change;
  Vertex tail, head, taild, headd=0, *id, *od;
  
  id=IN_DEG;
  od=OUT_DEG;
  lambda = INPUT_PARAM[0];
  oneexpl = 1.0-1.0/lambda;
  
  change = 0.0;
  FOR_EACH_TOGGLE(i) {
    echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 1;
    change += 4.0*echange;
    /* The above line may be an error -- should it be 2.0*echange? */
    taild = od[tail] + id[tail] + (echange - 1)/2;
    headd = od[head] + id[head] + (echange - 1)/2;
    if(taild!=0){
      change -= echange*(1.0-pow(oneexpl,(double)taild));
    }
    if(headd!=0){
      change -= echange*(1.0-pow(oneexpl,(double)headd));
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  CHANGE_STAT[0] = change;
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_gwb2
*****************/
D_CHANGESTAT_FN(d_gwb2){
  int i, echange=0;
  double alpha, oneexpa, change;
  Vertex tail, head, headd=0, *id;
  /* int nb1, nb2;

  nb2 = (int)INPUT_PARAM[0];
  nb1 = N_NODES - nb2; */

  id=IN_DEG;
  //od=OUT_DEG;
  change = 0.0;
  alpha = INPUT_PARAM[1];
  oneexpa = 1.0-exp(-alpha);
  
  FOR_EACH_TOGGLE(i) {
    echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 1;
    // taild = od[tail] + id[head] + (echange - 1)/2;
    headd = id[head] + (echange - 1)/2;
    // if(taild!=0){
      // change += echange*(1.0-pow(oneexpa,(double)taild));
    // }
    if(headd!=0){
      change += echange*(1.0-pow(oneexpa,(double)headd));
    }
    //Rprintf("tail %d head %d headd %d echange %d change %f\n", tail, head, headd, echange,change);
    // Rprintf(" tail %d head %d taild %d echange %d change %f\n", tail, head, taild, echange,change);
    // Rprintf(" od[tail] %d id[tail] %d od[head] %d id[head] %d\n", od[tail], id[tail], od[head], id[head]);
    TOGGLE_IF_MORE_TO_COME(i);
  }
  CHANGE_STAT[0] = change*exp(alpha);  
  // Rprintf("alpha  %f taild %d headd %d change %f\n", alpha, taild, headd, CHANGE_STAT[0]);
  UNDO_PREVIOUS_TOGGLES(i);
}

/********************  changestats:  H    ***********/
/*****************
 changestat: d_heideriandynamic
*****************/
D_CHANGESTAT_FN(d_heideriandynamic)  {
  long int nnodes;
  Vertex tail, head;
  int i, edgestate, edgestateheadtail;
  int edgestatep, edgestateheadtailp;
  
  nnodes = ((long int)(INPUT_PARAM[0]));
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i) {
    /*Get the initial edge state*/
    edgestate=IS_OUTEDGE(tail=TAIL(i), head=HEAD(i));
    edgestateheadtail = IS_OUTEDGE(head, tail);
    /*Get the prior edge state*/
    edgestatep  =(int)(INPUT_PARAM[1+(tail-1)+(head-1)*nnodes]);
    edgestateheadtailp=(int)(INPUT_PARAM[1+(head-1)+(tail-1)*nnodes]);

    // if(edgestatep != edgestateheadtailp){
      // if(edgestateheadtail!=edgestate){ CHANGE_STAT[0] += 1.0; }
    // }else{
      // if(edgestateheadtail==edgestate){ CHANGE_STAT[0] -= 1.0; }
    // }      
    if(edgestatep != edgestateheadtailp){
      if(edgestate != edgestateheadtail){
        CHANGE_STAT[0] += 1.0;
      }else{
        CHANGE_STAT[0] -= 1.0;
      }
    }
    // Rprintf("tail %d head %d edgestatep %d  edgestateheadtailp %d edgestate %d edgestateheadtail %d dstats %f\n",tail,head,edgestatep,edgestateheadtailp, edgestate, edgestateheadtail,  CHANGE_STAT[0] );
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
  // CHANGE_STAT[0] = 0.0;
  // Rprintf("nedges %d\n", N_EDGES );
  // Rprintf("dstats %f\n", CHANGE_STAT[0] );
}

/*****************
 changestat: d_hiertriad
*****************/
D_CHANGESTAT_FN(d_hiertriad) {
  Edge e;
  Vertex tail, head, node3;
  double pos3;
  int /* edgestate, */ i;
  
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i) {
    /* edgestate = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)); */
    tail = TAIL(i); head=HEAD(i);
    
    // Rprintf("tail %d head %d edgestate %d\n",tail,head, edgestate);
    STEP_THROUGH_OUTEDGES(head, e, node3) {
      if (IS_OUTEDGE(tail, node3)){
        pos3 = numposthree (node3, nwp);
        CHANGE_STAT[0] -= pos3;
        TOGGLE(tail, head);
        pos3 = numposthree (node3, nwp);
        CHANGE_STAT[0] += pos3;
        TOGGLE(tail, head);
      }
    }
    
    pos3 = numposthree (head, nwp);
    CHANGE_STAT[0] -= pos3;
    TOGGLE(tail, head);
    pos3 = numposthree (head, nwp);
    CHANGE_STAT[0] += pos3;
    TOGGLE(tail, head);
    // Rprintf("tail %d head %d pos3 %f ideg %d\n",tail,head, pos3, ideg);
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_hiertriaddegree
*****************/
D_CHANGESTAT_FN(d_hiertriaddegree) {
  Edge e;
  Vertex tail, head, node3, ideg;
  double pos3;
  int /* edgestate, */ i;
  
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i) {
    /* edgestate = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)); */
    tail = TAIL(i); head=HEAD(i);
    
    // Rprintf("tail %d head %d edgestate %d\n",tail,head, edgestate);
    STEP_THROUGH_OUTEDGES(head, e, node3) {
      if (IS_OUTEDGE(tail, node3)){
        ideg = IN_DEG[node3];
        pos3 = numposthree (node3, nwp);
        CHANGE_STAT[0] -= pos3*ideg;
        TOGGLE(tail, head);
        ideg = IN_DEG[node3];
        pos3 = numposthree (node3, nwp);
        CHANGE_STAT[0] += pos3*ideg;
        TOGGLE(tail, head);
      }
    }
    
    ideg = IN_DEG[head];
    pos3 = numposthree (head, nwp);
    CHANGE_STAT[0] -= pos3*ideg;
    TOGGLE(tail, head);
    ideg = IN_DEG[head];
    pos3 = numposthree (head, nwp);
    CHANGE_STAT[0] += pos3*ideg;
    TOGGLE(tail, head);    
    // Rprintf("tail %d head %d pos3 %f ideg %d\n",tail,head, pos3, ideg);
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
  numposthree:  called by d_hiertriad*
*****************/
double numposthree (Vertex head, Network *nwp) {
  Edge e, f;
  Vertex pos, node2, node3;
  double dpos;
  
  pos = 0;
  STEP_THROUGH_INEDGES(head, e, node2) {
    STEP_THROUGH_INEDGES(node2, f, node3) {
      if (IS_OUTEDGE(node3, head)) {++pos;} 
    }
    STEP_THROUGH_INEDGES(node2, f, node3) {
      if (IS_OUTEDGE(node3, head)) {++pos;} 
    }
  }
  dpos = pos / 2.0;
  return dpos;
}

/********************  changestats:  I    ***********/
/*****************
 changestat: d_icvar
*****************/
D_CHANGESTAT_FN(d_icvar)  {
  int i, edgestate, ichange, change;
  Vertex nnodes, tail, head, *id;
  
  id=IN_DEG;
  nnodes = N_NODES;
  
  change = 0;
  FOR_EACH_TOGGLE(i) {
      edgestate = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i));
      if(edgestate){
        ichange = -(2*(nnodes*(id[head]-1) - N_EDGES+1) + nnodes - 1);
      }else{
        ichange =   2*(nnodes* id[head]    - N_EDGES  ) + nnodes - 1;
      }
      // Rprintf("tail %d head %d nnodes %d  N_EDGES %d id[head] %d  ic %d\n",tail,head, nnodes,  N_EDGES, id[head], ichange);
      // change += edgestate ? (-ichange) : ichange;
      change += ichange;
      TOGGLE_IF_MORE_TO_COME(i);
  }
  CHANGE_STAT[0] = change*1.0/(nnodes*(nnodes-1));
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_idc
*****************/
D_CHANGESTAT_FN(d_idc)  {
  int i, edgestate, ichange, change;
  Vertex k, nnodes, maxidegree0, maxidegree1;
  Vertex tail, head, *id;
  
  id=IN_DEG;
  nnodes = N_NODES;
  
  change = 0;
  FOR_EACH_TOGGLE(i) {
      edgestate = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i));
      if(edgestate){
        maxidegree0 = id[head];
        maxidegree1 = id[head]-1;
        for (k=1; k<=nnodes; k++){
          if(id[k] > maxidegree0) {maxidegree0=id[k];}
          if(k != head && id[k] > maxidegree1) {maxidegree1=id[k];}
        }
        ichange = nnodes*(maxidegree1-maxidegree0) + 1;
      }else{
        maxidegree0 = 0;
        maxidegree1 = id[head]+1;
        for (k=1; k<=nnodes; k++){
          if(id[k] > maxidegree0) {maxidegree0=id[k];}
          if(id[k] > maxidegree1) {maxidegree1=id[k];}
        }
        ichange = nnodes*(maxidegree1-maxidegree0) - 1;
      }
      // Rprintf("tail %d head %d nnodes %d  N_EDGES %d id[head] %d  ic %d\n",tail,head, nnodes,  N_EDGES, id[head], ichange);
      change += ichange;
      TOGGLE_IF_MORE_TO_COME(i);
  }
  CHANGE_STAT[0] = change*1.0/((nnodes-1)*(nnodes-1));
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_intransitivedynamic
*****************/
D_CHANGESTAT_FN(d_intransitivedynamic)  {
  Edge e;
  long int nnodes;
  Vertex tail, head, node3;
  double change;
  int i, edgestate;
  
  nnodes = ((long int)(INPUT_PARAM[0]));
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i) {
    edgestate = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i));
    change = 0.0;
    STEP_THROUGH_OUTEDGES(head, e, node3) {
      if (node3 != tail){
        if (!IS_OUTEDGE(tail, node3)) {
          /* So form a tail -> head -> node3 intransitive triad */
          /* Check to see if a prior intransitive triad existed */
          if( INPUT_PARAM[1+(tail-1)+(head-1)*nnodes]==1 &&
            INPUT_PARAM[1+(head-1)+(node3-1)*nnodes]==1 &&
          INPUT_PARAM[1+(tail-1)+(node3-1)*nnodes]==0) {
            // Rprintf("tail %d head %d node3 %d nnodes %d\n",tail,head, node3, nnodes);
            change = change - 1.0;
          }
        }
      }
    }
    STEP_THROUGH_INEDGES(head, e, node3) {
      if (node3 != tail){
        if (IS_OUTEDGE(tail, node3)){
          /* So dissolve a tail -> node3 -> head intransitive triad */
          /*Check to see if a prior intransitive triad existed */
          if( INPUT_PARAM[1+(tail-1)+(head-1)*nnodes]==0 &&
            INPUT_PARAM[1+(tail-1)+(node3-1)*nnodes]==1 &&
          INPUT_PARAM[1+(node3-1)+(head-1)*nnodes]==1){
            change = change + 1.0;
          }
        }
      }
    }
    STEP_THROUGH_INEDGES(tail, e, node3) {
      if (node3 != head){
        if (IS_OUTEDGE(node3, head)) {
          /* So form a node3 -> tail -> head intransitive triad */
          /*Check to see if a prior intransitive triad existed */
          if( INPUT_PARAM[1+(tail-1)+(head-1)*nnodes]==1 &&
            INPUT_PARAM[1+(node3-1)+(tail-1)*nnodes]==1 &&
          INPUT_PARAM[1+(node3-1)+(head-1)*nnodes]==0) {
            change = change - 1.0;
          }
        }
      }
    }
    CHANGE_STAT[0] += edgestate ? -change : change;
    // Rprintf("tail %d head %d edgestate %d change %f\n",tail,head, edgestate, change);
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_intransitivity
*****************/
D_CHANGESTAT_FN(d_intransitivity) {
  int i, edgestate, a, b, c, d, e, edgecount, t300,
  t210, t120C, t120U, t120D, t201, t030C, t030T, t111U, 
  t111D, t021C, t021U, t021D, t102, t012 /*, t003 */ ;
  Vertex node3, tail, head;

  CHANGE_STAT[0] = 0.0;
  if (nwp->directed_flag) {
    // directed version
    FOR_EACH_TOGGLE(i) {
      edgestate = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i));
      t300 = 0;
      t210 = 0;
      t120C = 0;  t120U = 0;   t120D = 0;  t201 = 0;
      t030C = 0;  t030T = 0;   t111U = 0;  t111D = 0;
      t021C = 0;  t021U = 0;   t021D = 0;  t102 = 0;
      t012 = 0;
      
      if (MIN_OUTEDGE(head) != 0 || MIN_INEDGE(head) != 0 ||
      MIN_OUTEDGE(tail) != 0 || MIN_INEDGE(tail) != 0) {
        
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

              case 1:  { /* 021C, 021U, 021D, 102 */
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
              
              case 2: {   /* 030C, 030T, 111U, 111D */
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
              
              case 3: {  /* 120C, 120U, 120D, 201 */
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
              
              case 2:  { /* 030C, 030T, 111U, 111D */ 
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
                if (a == 1)
                {
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
                  } else
                  --t111U;
                }
              }
              break;
          
              case 4: {  /* 210 */
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
      
      // t003 = -(t300+t210+t120C+t120U+t120D+t201+t030C+t030T);
      // t003 = t003-(t111U+t111D+t021C+t021U+t021D+t102+t012);
      b = t021C+t030C+t111D+t111U+t120C+t201+t210;
      CHANGE_STAT[0] += edgestate ? -(double)b : (double)b;
      TOGGLE_IF_MORE_TO_COME(i);
    }
  }else{
    // undirected
    FOR_EACH_TOGGLE(i) {
      edgestate = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i));
      t300 = 0; t201 = 0; t102 = 0; t012 = 0;

      if (MIN_OUTEDGE(head) != 0 || MIN_INEDGE(head) != 0 ||
      MIN_OUTEDGE(tail) != 0 || MIN_INEDGE(tail) != 0 ) {
        
        /* ****** loop through node3 ****** */
        for (node3=1; node3 <= N_NODES; node3++) { 
          if (node3 != tail && node3 != head) {
            a = IS_OUTEDGE(MIN(node3,head), MAX(node3,head));
            b = IS_OUTEDGE(MIN(node3,tail), MAX(node3,tail));
            edgecount = (a + b);
            
            switch(edgecount) {
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
              
              case 2:  { /* 030C, 030T, 111U, 111D */
                ++t300;
                --t201;
              }
              break;
              
            }
          }  
        }    /* ******  move to next node3 ******** */
      } else 
      t102 = t102 + (N_NODES - 2);  
      
      /* t003 = -(t102+t201+t300); */
      b = t300;
      CHANGE_STAT[0] += edgestate ? -(double)b : (double)b;
      TOGGLE_IF_MORE_TO_COME(i);
    } // i loop
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/********************  changestats:  K    ***********/
/*****************
 changestat: d_kappa
*****************/
D_CHANGESTAT_FN(d_kappa)  {
  int i, j, echange=0;
  double nedges, change, ir0, fr0;
  Vertex tail, head, taild, headd=0, ik2, fk2, nnodes, *id, *od;
  
  id=IN_DEG;
  od=OUT_DEG;
  nnodes = N_NODES;
  
  change = 0.0;
  FOR_EACH_TOGGLE(i) {
    echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 1;
    ik2=0;
    for (j=1; j<=nnodes; j++) {      
      fk2 = od[j] + id[j];
      ik2 += fk2*(fk2-1);
    }
    taild = od[tail] + id[tail] + (echange-1)/2;
    headd = od[head] + id[head] + (echange-1)/2;
    fk2 = ik2 + echange*2*(taild+headd);
    nedges = (double)(N_EDGES);
    ir0 = (N_EDGES==0) ? 0.0 : (ik2*0.5/nedges);
    fr0 = (((N_EDGES)+echange)==0) ? 0.0 : (fk2*0.5/(nedges+echange));
    change += fr0 - ir0;
    // Rprintf("tail %d head %d nnodes %d nedges %f ik2 %d fk2 %d ir0 %f fr0 %f change %f\n",tail,head, nnodes,  nedges, ik2, fk2, ir0, fr0, change);
    TOGGLE_IF_MORE_TO_COME(i);
  }
  CHANGE_STAT[0] = change;
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_monopolymixmat
 ...has nothing to do with the board game; should be read 
 "mono, poly mixing matrix"
*****************/
D_CHANGESTAT_FN(d_monopolymixmat) {
  int edgestate, i;
  Edge e;
  Vertex tail, head, Fdeg, Mdeg, otherF, otherM;
  /* m/p means monogamous/polygamous; F/M means female/male */
  int mFmM, mFpM, pFmM;   
  /* NB: pFpM would be redundant since the total of all 4 is #edges */
  Vertex *od=OUT_DEG, *id=IN_DEG;

  CHANGE_STAT[0] = CHANGE_STAT[1] = CHANGE_STAT[2] = 0.0;
  FOR_EACH_TOGGLE(i) {
    edgestate = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i));
    Fdeg = od[tail];
    Mdeg = id[head];
    /* Calculate contribution from change of (F,M) edge only */
    mFmM = (Fdeg==0 && Mdeg==0) - (Fdeg==1 && Mdeg==1 && edgestate);
    mFpM = (Fdeg==0 && Mdeg>0) - (Fdeg==1 && Mdeg>1 && edgestate);
    pFmM = (Mdeg==0 && Fdeg>0) - (Mdeg==1 && Fdeg>1 && edgestate);
    /* Now calculate contribution from other partners of F or M */
    if(Fdeg - edgestate == 1) {/* Only case that concerns us */
      for(e = EdgetreeMinimum(nwp->outedges, tail);
      (otherM = nwp->outedges[e].value) != 0 && otherM == head;
      e = EdgetreeSuccessor(nwp->outedges, e)); /* This finds otherM */
      if (id[otherM] > 1) {
        mFpM += (Fdeg==1 ? -1 : 1);
      } else {
        mFmM += (Fdeg==1 ? -1 : 1);
        pFmM += (Fdeg==1 ? 1 : -1);
      }
    }
    if(Mdeg - edgestate == 1) {/* Similarly for Mdeg */
      for(e = EdgetreeMinimum(nwp->inedges, head);
      (otherF = nwp->inedges[e].value) != 0 && otherF == tail;
      e = EdgetreeSuccessor(nwp->inedges, e)); /* This finds otherF */
      if (od[otherF] > 1) { /* otherF is poly */
        pFmM += (Mdeg==1 ? -1 : 1);
      } else { /*otherF is mono */
        mFmM += (Mdeg==1 ? -1 : 1);
        mFpM += (Mdeg==1 ? 1 : -1);
      }
    }
    CHANGE_STAT[0] += (double) mFmM;
    CHANGE_STAT[1] += (double) mFpM;
    CHANGE_STAT[2] += (double) pFmM;    
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/********************  changestats:  L    ***********/

/********************  changestats:  M    ***********/

/********************  changestats:  N    ***********/

/********************  changestats:  O    ***********/

/********************  changestats:  R    ***********/

/********************  changestats:  S    ***********/
/*****************
 changestat: d_simmeliandynamic
*****************/
D_CHANGESTAT_FN(d_simmeliandynamic)  {
  Edge e;
  long int nnodes;
  Vertex tail, head, node3, change;
  int i, edgestate, edgestateheadtail;
  
  nnodes = ((long int)(INPUT_PARAM[0]));
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i) {
    /*Get the initial edge state*/
    edgestate=IS_OUTEDGE(tail=TAIL(i), head=HEAD(i));
    edgestateheadtail =  !IS_OUTEDGE(head, tail);
    
    if(!edgestateheadtail){
      /*Check to see if this will form a Simmelian */
      change = 0;
      STEP_THROUGH_OUTEDGES(head, e, node3) {
        if (IS_OUTEDGE(node3, tail) 
          && IS_OUTEDGE(tail, node3) 
        && IS_OUTEDGE(node3, head)){
          /*So this will form a tail,t,node3 Simmelian */
          /*Check to see if a prior Simmelian existed */
          if( INPUT_PARAM[1+(head-1)+(node3-1)*nnodes]==1 &&
            INPUT_PARAM[1+(tail-1)+(node3-1)*nnodes]==1 &&
          INPUT_PARAM[1+(node3-1)+(tail-1)*nnodes]==1 &&
          INPUT_PARAM[1+(head-1)+(tail-1)*nnodes]==1 &&
          INPUT_PARAM[1+(tail-1)+(head-1)*nnodes]==1 &&
          INPUT_PARAM[1+(node3-1)+(head-1)*nnodes]==1 ){
            ++change;
          }
        }
      }
      
      change = 6*change;
      CHANGE_STAT[0] += edgestate ? -(double)change : (double)change;
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_spatial
  Changescores for inhomogeneous Bernoulli graphs
*****************/
D_CHANGESTAT_FN(d_spatial) {
  int edgestate, i;
  Vertex tail, head, n;
  double llr,pb,alpha,gamma;

  n=N_NODES;
  pb=INPUT_PARAM[0];
  alpha=INPUT_PARAM[1];
  gamma=INPUT_PARAM[2];
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    edgestate = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i));
    /*Lower trianglization of adjacency matrix implies ith element corresponds
    to (row-1) + (col-1)*(n-1) - choose(col,2).*/
    llr = -log(((1+exp(pb))*pow(1+exp(alpha)*(INPUT_PARAM[2+(head-1)+(tail-1)*(n-1)-tail*(tail-1)/2]),exp(gamma)))/exp(pb)-1);
    CHANGE_STAT[0] += edgestate ? - llr : llr;
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/********************  changestats:  T    ***********/
/*****************
 changestat: d_transitivedynamic
*****************/
D_CHANGESTAT_FN(d_transitivedynamic)  {
  Edge e;
  long int nnodes;
  Vertex tail, head, node3;
  double change;
  int i, edgestate;
  
  nnodes = ((long int)(INPUT_PARAM[0]));
  CHANGE_STAT[0] = 0.0;
  FOR_EACH_TOGGLE(i) {
    edgestate = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i));
    change = 0.0;
    STEP_THROUGH_OUTEDGES(head, e, node3) {
      if (node3 != tail) {
        if (IS_OUTEDGE(tail, node3)) {
          /* So form a tail -> head -> node3 intransitive triad */
          /* Check to see if a prior intransitive triad existed */
          if( !(INPUT_PARAM[1+(tail-1)+(head-1)*nnodes]==1 &&
            INPUT_PARAM[1+(head-1)+(node3-1)*nnodes]==1 &&
          INPUT_PARAM[1+(tail-1)+(node3-1)*nnodes]==0)) {
            // Rprintf("tail %d head %d node3 %d nnodes %d\n",tail,head, node3, nnodes);
            change = change - 1.0;
          }
        }
      }
    }
    STEP_THROUGH_INEDGES(head, e, node3) {
      if (node3 != tail) {
        if (IS_OUTEDGE(tail, node3)) {
          /* So dissolve a tail -> node3 -> head intransitive triad */
          /*Check to see if a prior transitive triad existed */
          if( !(INPUT_PARAM[1+(tail-1)+(head-1)*nnodes]==0 &&
            INPUT_PARAM[1+(tail-1)+(node3-1)*nnodes]==1 &&
          INPUT_PARAM[1+(node3-1)+(head-1)*nnodes]==1)) {
            change = change + 1.0;
          }
        }
      }
    }
    STEP_THROUGH_INEDGES(tail, e, node3) {
      if (node3 != head){
        if (!IS_OUTEDGE(node3, head)){
          /* So form a node3 -> tail -> head intransitive triad */
          /*Check to see if a prior transitive triad existed */
          if( !(INPUT_PARAM[1+(tail-1)+(head-1)*nnodes]==1 &&
            INPUT_PARAM[1+(node3-1)+(tail-1)*nnodes]==1 &&
          INPUT_PARAM[1+(node3-1)+(head-1)*nnodes]==0)){
            change = change - 1.0;
          }
        }
      }
    }
    
    CHANGE_STAT[0] += edgestate ? -change : change;
    // Rprintf("tail %d head %d edgestate %d change %f\n",tail,head, edgestate, change);
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_transitivity
*****************/
D_CHANGESTAT_FN(d_transitivity) 
{
  int i, edgestate, a, b, c, d, e, edgecount, t300,
  t210, t120C, t120U, t120D, t201, t030C, t030T, t111U, 
  t111D, t021C, t021U, t021D, t102, t012, t003;
  Vertex node3, tail, head;

  CHANGE_STAT[0] = 0.0;

  if (nwp->directed_flag) {
    // directed version
    FOR_EACH_TOGGLE(i) {
      edgestate = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i));
      t300 = 0;
      t210 = 0;
      t120C = 0;  t120U = 0;   t120D = 0;  t201 = 0;
      t030C = 0;  t030T = 0;   t111U = 0;  t111D = 0;
      t021C = 0;  t021U = 0;   t021D = 0;  t102 = 0;
      t012 = 0;
      
      if (MIN_OUTEDGE(head) != 0 || MIN_INEDGE(head) != 0 || 
      MIN_OUTEDGE(tail) != 0 || MIN_INEDGE(tail) != 0) {
        
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
                 
                 case 3: {  /* 120C, 120U, 120D, 201 */
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

                 case 3:  { /* 201, 120D, 120U, 120C */
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
                     } else
                     --t111U;
                   }
                 }
                 break;
                 
                 case 4: {  /* 210 */
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
      } else 
      t012 = t012 + (N_NODES - 2);  
      
      t003 = -(t300+t210+t120C+t120U+t120D+t201+t030C+t030T);
      t003 = t003-(t111U+t111D+t021C+t021U+t021D+t102+t012);
      b = t003+t012+t021U+t021D+t102+t030T+t120U+t120D+t300;
      CHANGE_STAT[0] += edgestate ? -(double)b : (double)b;
      TOGGLE_IF_MORE_TO_COME(i);
    }
  } else {
    // undirected
    FOR_EACH_TOGGLE(i) {
      edgestate = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i));
      t300 = 0; t201 = 0; t102 = 0; t012 = 0;

      if (MIN_OUTEDGE(head) != 0 || MIN_INEDGE(head) != 0 ||
      MIN_OUTEDGE(tail) != 0 || MIN_INEDGE(tail) != 0 ) {
        
        /* ****** loop through node3 ****** */
        for (node3=1; node3 <= N_NODES; node3++) { 
          if (node3 != tail && node3 != head) {
            a = IS_OUTEDGE(MIN(node3,head), MAX(node3,head));
            b = IS_OUTEDGE(MIN(node3,tail), MAX(node3,tail));
            edgecount = (a + b);
            
            switch(edgecount)  {  
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

      t003 = (t102+t201+t300);
      b = -t300; 
      CHANGE_STAT[0] += edgestate ? -(double)b : (double)b;
      TOGGLE_IF_MORE_TO_COME(i);
    } // i loop
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_b1share
*****************/
D_CHANGESTAT_FN(d_b1share)  {
  Edge e, f;
  int i, j, echange;
  int L2tu;
  Vertex deg;
  Vertex tail, head, u, v;
  /* int nb2, nb1;

  nb1 = BIPARTITE;
  nb1 = (int)INPUT_PARAM[0];
  nb2 = (N_NODES) - nb1; */

  ZERO_ALL_CHANGESTATS(i);  
  FOR_EACH_TOGGLE(i) {
    echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 1;
    //
    // Next for b1 shared b2 counts
    STEP_THROUGH_INEDGES(head, e, u) {
      if (u != tail){
        L2tu=0;
        STEP_THROUGH_OUTEDGES(u, f, v) {
          if(IS_OUTEDGE(tail,v)) L2tu++;
  	    }
        for(j = 0; j < mtp->nstats; j++){
  	      deg = (Vertex)INPUT_PARAM[j+1];
          CHANGE_STAT[j] += ((L2tu + echange == deg)
          - (L2tu == deg));
  	    }
      }
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_b2share
*****************/
D_CHANGESTAT_FN(d_b2share)  {
  Edge e, f;
  int i, j, echange;
  int L2tu;
  Vertex deg;
  Vertex tail, head, u, v;
  /* int nb2, nb1;

  nb1 = (int)INPUT_PARAM[0];
  nb2 = (N_NODES) - nb1; */

  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
     echange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 1;
     //
     // Next for b2 shared b1 counts
     STEP_THROUGH_INEDGES(head, e, u) {
       if (u != tail){
         L2tu=0;
         STEP_THROUGH_OUTEDGES(u, v, f) {
           if(IS_OUTEDGE(tail,v)) L2tu++;
         }
         for(j = 0; j < mtp->nstats; j++){
           deg = (Vertex)INPUT_PARAM[j+1];
           CHANGE_STAT[j] += ((L2tu + echange == deg)
           - (L2tu == deg));
         }
       }
     }
     TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_gwb2share
****************/
D_CHANGESTAT_FN(d_gwb2share) {
  Edge e, f;
  int i, echange, ochange;
  int L2uh;
  Vertex tail, head, u, v;
  double alpha, oneexpa, cumchange;
  
  CHANGE_STAT[0] = 0.0;
  alpha = INPUT_PARAM[0];
  oneexpa = 1.0-exp(-alpha);
  
  FOR_EACH_TOGGLE(i) {
    cumchange=0.0;
    ochange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 0;
    echange = 2*ochange + 1;
    STEP_THROUGH_OUTEDGES(tail, e, u) {
      if (u != head){
        L2uh=ochange;
        STEP_THROUGH_INEDGES(u, v, f) {
          if(IS_OUTEDGE(MIN(v,head),MAX(v,head))) L2uh++;
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
 changestat: d_gwb1share
****************/
D_CHANGESTAT_FN(d_gwb1share) {
  Edge e, f;
  int i, echange, ochange;
  int L2tu;
  Vertex tail, head, u, v;
  double alpha, oneexpa, cumchange;
  
  CHANGE_STAT[0] = 0.0;
  alpha = INPUT_PARAM[0];
  oneexpa = 1.0-exp(-alpha);
  
  FOR_EACH_TOGGLE(i) {
    cumchange=0.0;
    ochange = IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) ? -1 : 0;
    echange = 2*ochange + 1;
    STEP_THROUGH_INEDGES(head, e, u) {
      if (u != tail){
        L2tu=ochange;
        STEP_THROUGH_OUTEDGES(u, v, f) {
          if(IS_OUTEDGE(MIN(v,tail),MAX(v,tail))) L2tu++;
        }
        cumchange += pow(oneexpa,(double)L2tu);
      }
    }
    cumchange  = echange*cumchange;
    (CHANGE_STAT[0]) += cumchange;
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

