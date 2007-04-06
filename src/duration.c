/*  This is a collection of functions used to calculate durations of
    edges and durations of overlaps.  They require a very specific set
    of input and return a very specific form of output.  DH 04-02-07 */

#include "duration.h"

/* These #defines are not really necessary but may make the code a bit
   easier to read.  They come at a price (the possibility of hard-to-track
   errors).   */
#define DMATRIX(a,b) (dmatrix[(a)+(offset)*(b)])
#define EDGE(a,b) (edge[(a)+(*nedge)*(b)])
#define CHANGE(a,b) (change[(a)+(*nchange)*(b)])
#define DISSOLVE(a,b) (dissolve[(a)+(*ndissolve)*(b)])
#define OMATRIX(a,b) (omatrix[(a)+(maxo)*(b)])


/**********************************************************/      

void DurationMatrix (int *nedge, int *edge, int *ntimestep, int *nfem, 
      int *ntotal, int *nchange, int *change, int *ndissolve, int *dissolve,
      int *dmatrix) {
  int row, j, k, f, m, time, offset = *nedge + *nchange;
  int i;

  /* Note:  This code assumes always that edges are listed in
     (female, male) order, where female < male.  */
  
  /* First, initialize dmatrix by putting in time-zero network */
  for (row=0; row<*nedge; row++) {
    AddNewDurationRow (dmatrix, row, EDGE(row,0), EDGE(row,1), 0, offset);
  } /* Note:  Value of i upon leaving loop is important */

  /* Next, step through time one click at a time */
  for (time=j=0; time<*ntimestep; time++) {
    for(; CHANGE(j,0) == time; j++) {
      f = CHANGE(j,1);
      m = CHANGE(j,2);
      for(k=row; !(DMATRIX(k, 0)==f && DMATRIX(k, 1)==m) && k>=0; k--);
      if (k>=0 && DMATRIX(k,3) == -1) { 
        /* We found a match for the (f, m) edge that must be ended */
        DMATRIX(k, 3) = time+1;
        DMATRIX(k, 4) = 1; /* non-Censored indicator */
        /* for(i = *ndissolve; DISSOLVE(i, 1)!=f && DISSOLVE(i,2)!=m && i>=0; i--);
        if (i<0) {
          Rprintf("Warning!  Dissolved edge (%d, %d) at time %d not contained in dissolve list\n",
          f,m,time);
        } */
      } else {
        AddNewDurationRow(dmatrix, row++, f, m, time+1, offset);
      }
    }
  }
}


void AddNewDurationRow (int *dmatrix, int row, int f, int m, int time, int offset) {
  DMATRIX(row, 0) = f;    /* female node number */
  DMATRIX(row, 1) = m;    /* male node number */
  DMATRIX(row, 2) = time; /* timestamp: edge begins */
  DMATRIX(row, 3) = -1;   /* timestamp:  edge ends or -1 if not yet known */
  DMATRIX(row, 4) = 0;    /* non-censoring indicator:  0=censored, 1=not */
}

void OverlapDurations (int *nedge, int *edge, int *ntimestep, int *nfem,
      int *ntotal, int *nchange, int *change, int *ndissolve, int *dissolve,
      int *maxoverlaps, int *omatrix) {
  Edge e, i, j, k, time, row, ne = *nedge, maxo = *maxoverlaps;
  Vertex f1, m1, f2, m2, nnodes, bipartite = *nfem;
  double *heads, *tails;
  TreeNode *ie, *oe;
  Network nw;

  /* Initialize by finding time-zero network overlaps */
  row=0;
  for(i=1; i < ne; i++) {
    f1=EDGE(i,0); 
    m1=EDGE(i,1);
    for (j=0; j < i; j++) {
      f2=EDGE(j,0); 
      m2=EDGE(j,1);
      if (f1 == f2 || m1 == m2) {
        AddNewOverlapRow(omatrix,row++,f1,m1,f2,m2,0,maxo);
      }              
    }
  }
  
  /* Set up a statnet Network using existing edgeTree code 
     Must coerce heads and tails to double, find nnodes */
  heads = (double *) malloc(sizeof(double) * ne);
  tails = (double *) malloc(sizeof(double) * ne);
  for (nnodes=1, i=0; i < ne; i++) {
    heads[i] = (double) EDGE(i, 0);
    m1 = EDGE(i, 1);
    if (m1 > nnodes) nnodes = m1; /* We assume m1 > f1 always */
    tails[i] = (double) m1;
  }
  nw = NetworkInitialize(heads, tails, ne, nnodes, 0, bipartite);
  free (heads);
  free (tails);
  ie=nw.inedges;
  oe=nw.outedges;

  /* Step through time one click at a time */
  for (time=j=0; time<*ntimestep; time++) {
    for(; CHANGE(j,0) == time; j++) {
      f1 = CHANGE(j,1);
      m1 = CHANGE(j,2);
      if (EdgetreeSearch(f1, m1, nw.outedges) != 0) { 
        /* Edge exists and must be removed. */
        /* Check whether there are any (zero-length) "concurrencies" with edges
        that begin AT THE CURRENT TIME but that haven't yet been encountered. */
        for (k=j+1; k < *nchange && CHANGE(k,0) == time; k++) {
          f2 = CHANGE(k, 1);
          m2 = CHANGE(k, 2);
          if ((f1==f2 || m1==m2) && EdgetreeSearch(f2,m2,nw.outedges)==0) {
            /* Any new concurrency begun this way will be immediately ended */
            AddNewOverlapRow(omatrix,row++,f1,m1,f1,m2,time+1,maxo); 
          } 
        }
        /* End all unfinished overlaps involving this edge */
        for(i=0; i <= row; i++) {
          if (OMATRIX(i, 5) == -1 &&
            ((f1 == OMATRIX(i, 0) && m1 == OMATRIX(i, 1)) ||
            (f1 == OMATRIX(i, 2) && m1 == OMATRIX(i, 3)))) {
              OMATRIX(i, 5) = time+1;  /* End time for overlap */
              OMATRIX(i, 6) = 1; /* non-censoring indicator */
          }
        }
      } else {  /* Edge does not exist and must be created. */
        /* First, add any new overlaps that are created. */
        /* Step through all the existing edges of f1 */
        /* (NB:  Since f<m, all edges of f1 are out-edges) */
        for(e = EdgetreeMinimum(oe, f1); 
        (m2 = oe[e].value) != 0;
        e = EdgetreeSuccessor(oe, e)) {
          AddNewOverlapRow(omatrix,row++,f1,m1,f1,m2,time+1,maxo);
        }
        /* Now step through all existing (in-)edges of m1 */
        for(e = EdgetreeMinimum(ie, m1);
        (f2 = ie[e].value) != 0;
        e = EdgetreeSuccessor(ie, e)) {
          AddNewOverlapRow(omatrix,row++,f1,m1,f2,m1,time+1,maxo);
        }
      }
      /* Finally, toggle (create or destroy) the edge */
      ToggleEdge(f1, m1, &nw); 
      if (row >= maxo) {
        Rprintf("Error! Value of maxoverlaps=%d too small\n",
                 maxo);
        return;
      }
    }
  }
  /* Free memory used by network object before returning */  
  NetworkDestroy (&nw);
}

void AddNewOverlapRow (int *omatrix, int row, int f1, int m1, 
int f2, int m2, int time, int maxo) {
  OMATRIX(row, 0) = f1; /* Female, edge 1 */
  OMATRIX(row, 1) = m1; /* Male, edge 1 */
  OMATRIX(row, 2) = f2; /* Female, edge 2 */
  OMATRIX(row, 3) = m2; /* Male, edge 2 */
  OMATRIX(row, 4) = time;  /* Start time for overlap */
  OMATRIX(row, 5) = -1; /* End time (-1 if not yet observed) */
  OMATRIX(row, 6) = 0;  /* noncensoring indicator: 0=censored */
}


/*****************
 void godfather_wrapper

 ...we'll make them an offer (of toggles) they can't refuse.
 This function takes a list of toggles, each with a time stamp,
 then produces a matrix of changestats (with one row for each unique
 time stamp value) that result from performing all the toggles at
 each time step.  For instance, one might use this function to 
 find the changestats that result from starting from an empty network
 and then adding all of the edges to make up an observed network of interest.
*****************/
void godfather_wrapper (double *heads, double *tails, double *dnedges,
                   double *dn, int *dflag, double *bipartite, 
                   int *nterms, char **funnames,
                   char **sonames, 
                   double *totalntoggles, double *timestamps, 
                   double *toggleheads, double *toggletails,
                   double *inputs, 
                   double *changestats, 
                   int *newnetwork, 
                   int *fVerbose, 
                   double *maxedges) {
  int directed_flag, maxtoggles, ntoggles;
  Vertex n_nodes, bip, *theads, *ttails, h, t;
  Edge e, i, j, n_edges, nmax, samplesize, nextedge, tnt;
  Network nw;
  Model *m;
  char *fn, *sn;
  ModelTerm *mtp;
  double *dstats, thistime; 
  
  n_nodes = (Vertex)*dn; /* coerce double *dn to type Vertex */
  n_edges = (Vertex)*dnedges; /* coerce double *dnedges to type Vertex */
  nmax = (Vertex)*maxedges; /* coerce double *maxedges to type Vertex */
  bip = (Vertex)*bipartite; /* coerce double *bipartite to type Vertex */    
  tnt = (Edge) *totalntoggles; /* coerce double *totalntoggles to type Edge */ 
  directed_flag = *dflag;

  for (i = 0; i < nmax; i++)
    newnetwork[i] = 0;
  m=ModelInitialize(*funnames, *sonames, inputs, *nterms);
  nw = NetworkInitialize(heads, tails, n_edges, n_nodes, directed_flag, bip);

  /* Figure out what is the maximum number of toggles at any time step and
  allocate an appropriate amount of memory. */
  thistime = timestamps[0];
  for(maxtoggles = samplesize = i = j = 1; i < tnt; i++) {
    if (timestamps[i] == thistime) {
      ++j; /* j counts how many timestamps match the current value. */
    } else { /* we have encountered a new timestamp value. */
      samplesize++;
      maxtoggles = maxtoggles < j? j : maxtoggles;
      j = 1;
      thistime = timestamps[i];
    }
  }
  theads = (Vertex *) malloc(sizeof(Vertex) * maxtoggles);
  ttails = (Vertex *) malloc(sizeof(Vertex) * maxtoggles);

  if (*fVerbose) {
    Rprintf("Total m->n_stats is %i; total samplesize is %d\n",
    m->n_stats,samplesize);
    Rprintf("maxedges = %ld, maxtoggles = %ld, samplesize = %ld, totalntoggles = %ld\n",
    nmax, maxtoggles, samplesize, tnt);
  }
  
  /*********************
  changestats are modified in groups of m->n_stats, and they
  reflect the CHANGE in the values of the statistics from the
  original (observed) network.  Thus, when we begin, the initial 
  values of the first group of m->n_stats changestats should 
  all be zero
  *********************/
  for (j=0; j < m->n_stats; j++)
    changestats[j] = 0.0;
  mtp = m->termarray;
  
  /* Now start obtaining change statistics */
  thistime = timestamps[0];
  theads[0] = toggleheads[0];
  ttails[0] = toggletails[0];
  for(i = ntoggles = 1; i <= tnt; i++) {
    if (i<tnt && timestamps[i] == thistime) {
      theads[ntoggles] = toggleheads[i];
      ttails[ntoggles] = toggletails[i];
      ++ntoggles; /* ntoggles counts how many timestamps match the current value. */
    } else { /* we have encountered a new timestamp value. */
      mtp = m->termarray;
      dstats = m->workspace;
      /* Time to calculate the change statistics for these toggles */
      /* Set current vector of stats equal to previous vector */
      for (j=0; j<m->n_stats; j++){
        changestats[j+m->n_stats] = changestats[j];
      }
      changestats += m->n_stats;
      /* This then adds the change statistics to these values */
      dstats = changestats;
      for (j=0; j < m->n_terms; j++) {  /* loop through each term */
        mtp->dstats = dstats;
        (*(mtp->func))(ntoggles, theads, ttails, mtp, &nw);  /* Call d_xxx function */
        dstats += (mtp++)->nstats;
      }
      /* Make proposed toggles for real this time) */
      for (j=0; j < ntoggles; j++){
        ToggleEdge(theads[j], ttails[j], &nw);
      }
      /* Finished with this timestamp; go on to next one */
      if (i<tnt) {
        ntoggles = 1;
        thistime = timestamps[i];
        theads[0] = toggleheads[i];
        ttails[0] = toggletails[i];
      }
    }
  }
  /* record new generated network to pass back to R */
  for (h=nextedge=1; h<=n_nodes; h++) {
    for(e = EdgetreeMinimum(nw.outedges, h);
    (t=nw.outedges[e].value) != 0 && nextedge < nmax;
    e = EdgetreeSuccessor(nw.outedges, e)) {
      newnetwork[nextedge++] = h;
      newnetwork[nextedge++] = t;
    }
  }
  if (nextedge>=nmax) { 
    Rprintf("Error!  The value of maxedges was not set high enough in ergm.godfather\n");
  }
  newnetwork[0]=nextedge;

  /* Clean up and return */
  free(theads);
  free(ttails);
  ModelDestroy(m);
  NetworkDestroy(&nw);
}

/*****************
 void mixingmatrix_hack

 This is a modified version of godfather_wrapper.  It has a very specific
 job and this is not a good way to go about this, but there does not appear
 to be any appropriate change statistic in place.
 The needed change statistics are the entries in a 2x2 mixing matrix that
 categorizes all of the edges in the network into one of four categories: 
 edges between two monogamous people, edges between a monogamous female and
 a polygamous male, edges between a monogamous male and a polygamous female,
 and edges between two polygamous people.

 THIS IS A TEMPORARY SOLUTION!
 
*****************/
void mixingmatrix_hack (double *heads, double *tails, double *dnedges,
                   double *dn, int *dflag, double *bipartite, 
                   int *nterms, char **funnames,
                   char **sonames, 
                   double *totalntoggles, double *timestamps, 
                   double *toggleheads, double *toggletails,
                   double *inputs, 
                   double *changestats, 
                   int *newnetwork, 
                   int *fVerbose, 
                   double *maxedges) {
  int directed_flag, maxtoggles, ntoggles;
  Vertex n_nodes, bip, *theads, *ttails, h, t;
  Edge e, i, j, n_edges, nmax, samplesize, nextedge, tnt;
  Network nw;
  Model *m;
  char *fn, *sn;
  ModelTerm *mtp;
  double *dstats, thistime; 
  
  n_nodes = (Vertex)*dn; /* coerce double *dn to type Vertex */
  n_edges = (Vertex)*dnedges; /* coerce double *dnedges to type Vertex */
  nmax = (Vertex)*maxedges; /* coerce double *maxedges to type Vertex */
  bip = (Vertex)*bipartite; /* coerce double *bipartite to type Vertex */    
  tnt = (Edge) *totalntoggles; /* coerce double *totalntoggles to type Edge */ 
  directed_flag = *dflag;

  for (i = 0; i < nmax; i++)
    newnetwork[i] = 0;
  m=ModelInitialize(*funnames, *sonames, inputs, *nterms);
  nw = NetworkInitialize(heads, tails, n_edges, n_nodes, directed_flag, bip);

  /* Figure out what is the maximum number of toggles at any time step and
  allocate an appropriate amount of memory. */
  thistime = timestamps[0];
  for(maxtoggles = samplesize = i = j = 1; i < tnt; i++) {
    if (timestamps[i] == thistime) {
      ++j; /* j counts how many timestamps match the current value. */
    } else { /* we have encountered a new timestamp value. */
      samplesize++;
      maxtoggles = maxtoggles < j? j : maxtoggles;
      j = 1;
      thistime = timestamps[i];
    }
  }
  theads = (Vertex *) malloc(sizeof(Vertex) * maxtoggles);
  ttails = (Vertex *) malloc(sizeof(Vertex) * maxtoggles);

  if (*fVerbose) {
    Rprintf("Total m->n_stats is %i; total samplesize is %d\n",
    m->n_stats,samplesize);
    Rprintf("maxedges = %ld, maxtoggles = %ld, samplesize = %ld, totalntoggles = %ld\n",
    nmax, maxtoggles, samplesize, tnt);
  }
  
  /*********************
  changestats are modified in groups of m->n_stats, and they
  reflect the CHANGE in the values of the statistics from the
  original (observed) network.  Thus, when we begin, the initial 
  values of the first group of m->n_stats changestats should 
  all be zero
  *********************/
  for (j=0; j < m->n_stats; j++)
    changestats[j] = 0.0;
  mtp = m->termarray;
  
  /* Now start obtaining change statistics */
  thistime = timestamps[0];
  theads[0] = toggleheads[0];
  ttails[0] = toggletails[0];
  for(i = ntoggles = 1; i <= tnt; i++) {
    if (i<tnt && timestamps[i] == thistime) {
      theads[ntoggles] = toggleheads[i];
      ttails[ntoggles] = toggletails[i];
      ++ntoggles; /* ntoggles counts how many timestamps match the current value. */
    } else { /* we have encountered a new timestamp value. */
      mtp = m->termarray;
      dstats = m->workspace;
      /* Time to calculate the change statistics for these toggles */
      /* Set current vector of stats equal to previous vector */
      for (j=0; j<m->n_stats; j++){
        changestats[j+m->n_stats] = changestats[j];
      }
      changestats += m->n_stats;
      /* This then adds the change statistics to these values */
      dstats = changestats;
      for (j=0; j < m->n_terms; j++) {  /* loop through each term */
        mtp->dstats = dstats;
        dstats[0] = 
        dstats[1] =         
        dstats[2] =         
        dstats[3] =         
        (*(mtp->func))(ntoggles, theads, ttails, mtp, &nw);  /* Call d_xxx function */
        dstats += (mtp++)->nstats;
      }
      /* Make proposed toggles for real this time) */
      for (j=0; j < ntoggles; j++){
        ToggleEdge(theads[j], ttails[j], &nw);
      }
      /* Finished with this timestamp; go on to next one */
      if (i<tnt) {
        ntoggles = 1;
        thistime = timestamps[i];
        theads[0] = toggleheads[i];
        ttails[0] = toggletails[i];
      }
    }
  }
  /* record new generated network to pass back to R */
  for (h=nextedge=1; h<=n_nodes; h++) {
    for(e = EdgetreeMinimum(nw.outedges, h);
    (t=nw.outedges[e].value) != 0 && nextedge < nmax;
    e = EdgetreeSuccessor(nw.outedges, e)) {
      newnetwork[nextedge++] = h;
      newnetwork[nextedge++] = t;
    }
  }
  if (nextedge>=nmax) { 
    Rprintf("Error!  The value of maxedges was not set high enough in ergm.godfather\n");
  }
  newnetwork[0]=nextedge;

  /* Clean up and return */
  free(theads);
  free(ttails);
  ModelDestroy(m);
  NetworkDestroy(&nw);
}




