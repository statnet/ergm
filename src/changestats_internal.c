/***************************************************************
 changestats internal
***************************************************************/

/*****************
 changestat: d_b1degree_edgecov
*****************/

CHANGESTAT_FN(d_b1degree_edgecov) { 
  int i, j, k, echange, n1, n2, edgecovval, min, max, mid, nedges;
  Vertex h1, t1, h2, t2, b1deg, d;

  n1 = BIPARTITE;
  n2 = N_NODES - BIPARTITE;
  nedges = (N_INPUT_PARAMS - N_CHANGE_STATS)/2;
  for (i=0; i < N_CHANGE_STATS; i++) 
    CHANGE_STAT[i] = 0.0;  
  FOR_EACH_TOGGLE(i) {
    b1deg = 0;
    h1 = heads[i];
    t1 = tails[i];
    echange = IS_OUTEDGE(h1, t1) ? -1 : 1;
    edgecovval = 0;
    min=N_CHANGE_STATS;          /*  Begin: Determine whether this tie exists in edgecov*/
    max=nedges+min-1;
    while (max >= min) {
      mid = (min + max)/2;
      h2 = INPUT_PARAM[mid+N_CHANGE_STATS-1];
      t2 = INPUT_PARAM[mid+nedges+N_CHANGE_STATS-1];
      if (h1<h2 || ((h1==h2)&&(t1<t2))) { /* Move search window down */
        max = mid-1;
      }
      else if (h1>h2 || ((h1==h2)&&(t1>t2))) { /* Move search window up */
        min = mid+1;
      }
      else {
        edgecovval = 1;
        break;
      }
    }                            /*  End: Determine whether this tie exists in edgecov*/
    if (edgecovval == 1) {       /*  Begin: determine current edgecov-degree of node of interest*/
      for (k=n1+1; k <= N_NODES; k++) 
        if (IS_OUTEDGE(h1,k)) {
          min=N_CHANGE_STATS;
          max=nedges+min-1;
          while (max >= min) {
            mid = (min + max)/2;
            h2 = INPUT_PARAM[mid+N_CHANGE_STATS-1];
            t2 = INPUT_PARAM[mid+nedges+N_CHANGE_STATS-1];
            if (h1<h2 || ((h1==h2)&&(k<t2))) { /* Move search window down */
              max = mid-1;
            }
            else if (h1>h2 || ((h1==h2)&&(k>t2))) { /* Move search window up */
              min = mid+1;
            }
            else {
              b1deg ++;
              break;
            }
          }
        }                      
    }                            /*  End: determine current edgecov-degree of node of interest*/
    for(j = 0; j < N_CHANGE_STATS; j++) {
      d = (Vertex)(INPUT_PARAM[j]);
      CHANGE_STAT[j] += (b1deg + echange*edgecovval == d) - (b1deg == d);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_b2degree_edgecov
*****************/

CHANGESTAT_FN(d_b2degree_edgecov) { 
  int i, j, k, echange, n1, n2, edgecovval, min, max, mid, nedges;
  Vertex h1, t1, h2, t2, b2deg, d;

  n1 = BIPARTITE;
  n2 = N_NODES - BIPARTITE;
  nedges = (N_INPUT_PARAMS - N_CHANGE_STATS)/2;  
  for (i=0; i < N_CHANGE_STATS; i++)
    CHANGE_STAT[i] = 0.0;  
  FOR_EACH_TOGGLE(i) {
    b2deg = 0;
    h1 = heads[i];
    t1 = tails[i];
    echange = IS_OUTEDGE(h1, t1) ? -1 : 1;
    edgecovval = 0;
    min=N_CHANGE_STATS;          /*  Begin: Determine whether this tie exists in edgecov*/
    max=nedges+min-1;
    while (max >= min) {
      mid = (min + max)/2;
      h2 = INPUT_PARAM[mid+N_CHANGE_STATS-1];
      t2 = INPUT_PARAM[mid+nedges+N_CHANGE_STATS-1];
      if (h1<h2 || ((h1==h2)&&(t1<t2))) { /* Move search window down */
        max = mid-1;
      }
      else if (h1>h2 || ((h1==h2)&&(t1>t2))) { /* Move search window up */
        min = mid+1;
      }
      else {
        edgecovval = 1;
        break;                  
      }
    }                            /*  End: Determine whether this tie exists in edgecov*/
    if (edgecovval == 1) {       /*  Begin: determine current edgecov-degree of node of interest*/
      for (k=1; k <= n1; k++) 
        if (IS_OUTEDGE(k,t1)) {
          min=N_CHANGE_STATS;
          max=nedges+min-1;
          while (max >= min) {
            mid = (min + max)/2;
            h2 = INPUT_PARAM[mid+N_CHANGE_STATS-1];
            t2 = INPUT_PARAM[mid+nedges+N_CHANGE_STATS-1];
            if (k<h2 || ((k==h2)&&(t1<t2))) { /* Move search window down */
              max = mid-1;
            }
            else if (k>h2 || ((k==h2)&&(t1>t2))) { /* Move search window up */
              min = mid+1;
            }
            else {
              b2deg ++;
              break;
            }
          }
        }                      
    }                            /*  End: determine current edgecov-degree of node of interest*/
    for(j = 0; j < N_CHANGE_STATS; j++) {
      d = (Vertex)(INPUT_PARAM[j]);
      CHANGE_STAT[j] += (b2deg + echange*edgecovval == d) - (b2deg == d);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}


/*****************
 changestat: d_b1mindegree
*****************/
CHANGESTAT_FN(d_b1mindegree) { 
  int i, j, echange;
  Vertex b1, b1deg, d;

  for (i=0; i < N_CHANGE_STATS; i++) 
    CHANGE_STAT[i] = 0.0;  
  FOR_EACH_TOGGLE(i) {
    b1 = heads[i];
    echange = IS_OUTEDGE(b1, tails[i]) ? -1 : 1;
    b1deg = OUT_DEG[b1];
    for(j = 0; j < N_CHANGE_STATS; j++) {
      d = (Vertex)(INPUT_PARAM[j]);
      CHANGE_STAT[j] += (b1deg + echange >= d) - (b1deg >= d);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}

/*****************
 changestat: d_b2mindegree
*****************/
CHANGESTAT_FN(d_b2mindegree) { 
  /* It is assumed that in this bipartite network, the only edges are
  of the form (b1, b2), where b1 is always strictly less
  than b2.  In other words, the degree of a b1 is equivalent
  to its outdegree and the degree of a b2 is equivalent to its
  indegree.
  */
  int i, j, echange;
  Vertex b1, b2, b2deg, d, *id;

  id=IN_DEG;
  for (i=0; i < N_CHANGE_STATS; i++) 
    CHANGE_STAT[i] = 0.0;  
  for (i=0; i<ntoggles; i++) {      
    echange=(EdgetreeSearch(b1=heads[i], b2=tails[i], nwp->outedges)==0) ? 1 : -1;
    b2deg = id[b2];
    for(j = 0; j < N_CHANGE_STATS; j++) {
      d = (Vertex)(INPUT_PARAM[j]);
      CHANGE_STAT[j] += (b2deg + echange >= d) - (b2deg >= d);
    }
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}


/*****************
 changestat: d_b1mindegree_edgecov
*****************/
CHANGESTAT_FN(d_b1mindegree_edgecov) { 
  int i, j, k, echange, n1, n2, edgecovval, min, max, mid, nedges;
  Vertex h1, t1, h2, t2, b1deg, d;

  n1 = BIPARTITE;
  n2 = N_NODES - BIPARTITE;
  nedges = (N_INPUT_PARAMS - N_CHANGE_STATS)/2;
  for (i=0; i < N_CHANGE_STATS; i++) 
    CHANGE_STAT[i] = 0.0;  
  FOR_EACH_TOGGLE(i) {
    b1deg = 0;
    h1 = heads[i];
    t1 = tails[i];
    echange = IS_OUTEDGE(h1, t1) ? -1 : 1;
    edgecovval = 0;
    min=N_CHANGE_STATS;          /*  Begin: Determine whether this tie exists in edgecov*/
    max=nedges+min-1;
    while (max >= min) {
      mid = (min + max)/2;
      h2 = INPUT_PARAM[mid+N_CHANGE_STATS-1];
      t2 = INPUT_PARAM[mid+nedges+N_CHANGE_STATS-1];
      if (h1<h2 || ((h1==h2)&&(t1<t2))) { /* Move search window down */
        max = mid-1;
      }
      else if (h1>h2 || ((h1==h2)&&(t1>t2))) { /* Move search window up */
        min = mid+1;
      }
      else {
        edgecovval = 1;
        break;
      }
    }                            /*  End: Determine whether this tie exists in edgecov*/
    if (edgecovval == 1) {       /*  Begin: determine current edgecov-degree of node of interest*/
      for (k=n1+1; k <= N_NODES; k++) 
        if (IS_OUTEDGE(h1,k)) {
          min=N_CHANGE_STATS;
          max=nedges+min-1;
          while (max >= min) {
            mid = (min + max)/2;
            h2 = INPUT_PARAM[mid+N_CHANGE_STATS-1];
            t2 = INPUT_PARAM[mid+nedges+N_CHANGE_STATS-1];
            if (h1<h2 || ((h1==h2)&&(k<t2))) { /* Move search window down */
              max = mid-1;
            }
            else if (h1>h2 || ((h1==h2)&&(k>t2))) { /* Move search window up */
              min = mid+1;
            }
            else {
              b1deg ++;
              break;
            }
          }
        }                      
    }                            /*  End: determine current edgecov-degree of node of interest*/
    for(j = 0; j < N_CHANGE_STATS; j++) {
      d = (Vertex)(INPUT_PARAM[j]);
      CHANGE_STAT[j] += (b1deg + echange*edgecovval >= d) - (b1deg >= d);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}


/*****************
 changestat: d_b2mindegree_edgecov
*****************/

CHANGESTAT_FN(d_b2mindegree_edgecov) { 
  int i, j, k, echange, n1, n2, edgecovval, min, max, mid, nedges;
  Vertex h1, t1, h2, t2, b2deg, d;

  n1 = BIPARTITE;
  n2 = N_NODES - BIPARTITE;
  nedges = (N_INPUT_PARAMS - N_CHANGE_STATS)/2;  
  for (i=0; i < N_CHANGE_STATS; i++)
    CHANGE_STAT[i] = 0.0;  
  FOR_EACH_TOGGLE(i) {
    b2deg = 0;
    h1 = heads[i];
    t1 = tails[i];
    echange = IS_OUTEDGE(h1, t1) ? -1 : 1;
    edgecovval = 0;
    min=N_CHANGE_STATS;          /*  Begin: Determine whether this tie exists in edgecov*/
    max=nedges+min-1;
    while (max >= min) {
      mid = (min + max)/2;
      h2 = INPUT_PARAM[mid+N_CHANGE_STATS-1];
      t2 = INPUT_PARAM[mid+nedges+N_CHANGE_STATS-1];
      if (h1<h2 || ((h1==h2)&&(t1<t2))) { /* Move search window down */
        max = mid-1;
      }
      else if (h1>h2 || ((h1==h2)&&(t1>t2))) { /* Move search window up */
        min = mid+1;
      }
      else {
        edgecovval = 1;
        break;                  
      }
    }                            /*  End: Determine whether this tie exists in edgecov*/
    if (edgecovval == 1) {       /*  Begin: determine current edgecov-degree of node of interest*/
      for (k=1; k <= n1; k++) 
        if (IS_OUTEDGE(k,t1)) {
          min=N_CHANGE_STATS;
          max=nedges+min-1;
          while (max >= min) {
            mid = (min + max)/2;
            h2 = INPUT_PARAM[mid+N_CHANGE_STATS-1];
            t2 = INPUT_PARAM[mid+nedges+N_CHANGE_STATS-1];
            if (k<h2 || ((k==h2)&&(t1<t2))) { /* Move search window down */
              max = mid-1;
            }
            else if (k>h2 || ((k==h2)&&(t1>t2))) { /* Move search window up */
              min = mid+1;
            }
            else {
              b2deg ++;
              break;
            }
          }
        }                      
    }                            /*  End: determine current edgecov-degree of node of interest*/
    for(j = 0; j < N_CHANGE_STATS; j++) {
      d = (Vertex)(INPUT_PARAM[j]);
      CHANGE_STAT[j] += (b2deg + echange*edgecovval >= d) - (b2deg >= d);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}


