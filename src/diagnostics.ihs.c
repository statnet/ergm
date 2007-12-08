/*  This is a collection of functions used to calculate diagnostic 
    statistics for dynamic networks. */

#include "diagnostics.ihs.h"

/* These #defines are not really necessary but may make the code a bit
   easier to read.  They come at a price (the possibility of hard-to-track
   errors).   */
#define DMATRIX(a,b) (dmatrix[(a)+(offset)*(b)])
#define EDGE(a,b) (edge[(a)+(*nedge)*(b)])
#define CHANGE(a,b) (change[(a)+(*nchange)*(b)])
#define OMATRIX(a,b) (omatrix[(a)+(maxo)*(b)])
#define DEGMIXMAT(a,b) (degmixmat[(a)+(*nnodes)*(b)])


/* Helper functions defined inline. */

R_INLINE void AddNewDurationRow (int *dmatrix, int row, int h, int t, int time, int offset) {
  DMATRIX(row, 0) = h;    /* head node number */
  DMATRIX(row, 1) = t;    /* tail node number */
  DMATRIX(row, 2) = time; /* timestamp: edge begins */
  DMATRIX(row, 3) = -1;   /* timestamp:  edge ends or -1 if not yet known */
  DMATRIX(row, 4) = 0;    /* non-censoring indicator:  0=censored, 1=not */
}

R_INLINE void AddNewOverlapRow (int *omatrix, int row, int f1, int m1, 
int f2, int m2, int time1, int time2, int maxo) {
  OMATRIX(row, 0) = f1; /* Female, edge 1 */
  OMATRIX(row, 1) = m1; /* Male, edge 1 */
  OMATRIX(row, 2) = f2; /* Female, edge 2 */
  OMATRIX(row, 3) = m2; /* Male, edge 2 */
  OMATRIX(row, 4) = time1;  /* Start time for edge 1 */
  OMATRIX(row, 5) = time2;  /* Start time for edge 2 */
  OMATRIX(row, 6) = -1; /* End time (-1 if not yet observed) */
  OMATRIX(row, 7) = 0;  /* noncensoring indicator: 0=censored, nonzero value
                           gives edge (1 or 2) that ended first (or 3 if they
                           ended simultaneously). */
}

/**********************************************************/      


void DurationMatrix (int *nedge, int *edge, int *ntimestep,
      int *ntotal, int *nchange, int *change,
      int *dmatrix) {
  int row, j, k, h, t, time, offset = *nedge + *nchange;

  /* Note:  This code assumes always that edges are listed in
     (head, tail) order, where, for bipartite and underected networks, head < tail.  */
  
  /* First, initialize dmatrix by putting in time-zero network */
  for (row=0; row<*nedge; row++) {
    AddNewDurationRow (dmatrix, row, EDGE(row,0), EDGE(row,1), 0, offset);
  } /* Note:  Value of i upon leaving loop is important */

  /* Next, step through time one click at a time */
  for (time=1,j=0; time<=*ntimestep; time++) {
    for(; CHANGE(j,0) == time; j++) {
      h = CHANGE(j,1);
      t = CHANGE(j,2);
      for(k=row; !(DMATRIX(k, 0)==h && DMATRIX(k, 1)==t) && k>=0; k--);
      if (k>=0 && DMATRIX(k,3) == -1) { 
        /* We found a match for the (h, t) edge that must be ended */
        DMATRIX(k, 3) = time;
        DMATRIX(k, 4) = 1; /* non-Censored indicator */
        /* int i;
        for(i = *ndissolve; DISSOLVE(i, 1)!=h && DISSOLVE(i,2)!=t && i>=0; i--);
        if (i<0) {
          Rprintf("Warning!  Dissolved edge (%d, %d) at time %d not contained in dissolve list\n",
          h,t,time);
        } */
      } else {
        AddNewDurationRow(dmatrix, row++, h, t, time, offset);
      }
    }
  }
}

void OverlapDurations (int *nnodes,
      int *nedge, int *edge, int *ntimestep, int *nfem,
      int *ntotal, int *nchange, int *change,
      int *maxoverlaps, int *omatrix) {
  Edge e, i, j, k, row, ne = *nedge;
  int f1, m1, time, f2, m2, match;
  int bipartite = *nfem, maxo=*maxoverlaps;
  double t1;
  WtTreeNode *ie, *oe;
  WtNetwork nw;

  /* Initialize by finding time-zero network overlaps */
  row=0;
  for(i=1; i < ne; i++) {
    f1=EDGE(i,0); 
    m1=EDGE(i,1);
    for (j=0; j < i; j++) {
      f2=EDGE(j,0); 
      m2=EDGE(j,1);
      if (f1 == f2 || m1 == m2) {
        AddNewOverlapRow(omatrix,row++,f1,m1,f2,m2,0,0,maxo);
      }              
    }
  }
  double *starttimes = (double *) malloc(sizeof(double) * ne);
  for (i=0; i < ne; i++) {
    starttimes[i] = 0.0;
  }

  /* R's serialization of matrixes is column-major, so this works: */
  nw = WtNetworkInitialize(edge, edge+*nedge, starttimes, ne, *nnodes, 0, bipartite);
  free (starttimes);
  ie=nw.inedges;
  oe=nw.outedges;

  /* Step through time one click at a time */
  for (time=1,j=0; time<=*ntimestep; time++) {
    for(; CHANGE(j,0) == time; j++) {
      t1 = (double) CHANGE(j,0);        
      f1 = CHANGE(j,1);
      m1 = CHANGE(j,2);
      if ((e=WtEdgetreeSearch(f1, m1, nw.outedges)) != 0) { 
        /* Edge exists and must be removed. */
        /* Check whether there are any (zero-length) "concurrencies" with edges
        that begin AT THE CURRENT TIME but that haven't yet been encountered. */
        for (k=j+1; k < *nchange && CHANGE(k,0) == time; k++) {
          f2 = CHANGE(k, 1);
          m2 = CHANGE(k, 2);
          if ((f1==f2||m1==m2) && WtEdgetreeSearch(f2,m2,nw.outedges)==0) {
            /* Edge1 is ending, edge2 beginning; create a 0-length conc. */
            /* Any new concurrency begun this way will be immediately ended */
            t1 = nw.outedges[e].weight;
            AddNewOverlapRow(omatrix, row++, f1, m1, f2, m2,
                             (int) t1, time, maxo);
          }
        }
        /* End all unfinished overlaps involving this edge */
        for(i=0; i <= row; i++) {
          match = (f1==OMATRIX(i, 0) && m1 == OMATRIX(i,1))? 1 : 0;
          match += (f1==OMATRIX(i, 2) && m1 == OMATRIX(i,3)) ? 2 : 0;
          if (match>0 && OMATRIX(i,6) == -1) { /* overlap has ended */
            OMATRIX(i,6) = time; /* End time */
            OMATRIX(i,7) = match; /* Which edge has ended? (1 or 2) */
            /* Check to see if the other edge ends at the 
            current time.  If so, then set OMATRIX(i,7)=3 */
            f2 = OMATRIX(i, 4-2*match); /* either OMATRIX(i, 2) or (i, 0) */
            m2 = OMATRIX(i, 5-2*match); /* either OMATRIX(i, 3) or (i, 1) */
            for (k=j+1; k < *nchange && CHANGE(k,0) == time; k++) {
              if (CHANGE(k,1) == f2 && CHANGE(k,2) == m1) { 
                OMATRIX(i,7)=3; /* Both edges ended simultaneously */
              }
            }
          }
        }
      } else {  /* Edge does not exist and must be created. */
        /* First, add any new overlaps that are created. */
        /* Step through all the existing edges of f1 */
        /* (NB:  Since f<m, all edges of f1 are out-edges) */
        for(e = WtEdgetreeMinimum(oe, f1); 
        (m2 = oe[e].value) != 0;
        e = WtEdgetreeSuccessor(oe, e)) {
          AddNewOverlapRow(omatrix, row++, f1, m1, f1, m2, 
                           time, (int)oe[e].weight, maxo);
        }
        /* Now step through all existing (in-)edges of m1 */
        for(e = WtEdgetreeMinimum(ie, m1);
        (f2 = ie[e].value) != 0;
        e = WtEdgetreeSuccessor(ie, e)) {
          AddNewOverlapRow(omatrix,row++,f1,m1,f2,m1,
                           time, (int) ie[e].weight, maxo);
        }
      }
      /* Finally, toggle (create or destroy) the edge */
        WtToggleEdge(f1, m1, 1.0+t1, &nw); 
      if (row >= maxo) {
        Rprintf("Error! Value of maxoverlaps=%d too small\n",
                 maxo);
        return;
      }
    }
  }
  /* Free memory used by network object before returning */  
  WtNetworkDestroy(&nw);
}

void DegreeMixMatrix (int *nnodes,
      int *nedge, int *edge, int *ntimestep, int *nfem,
      int *ntotal, int *nchange, int *change,
      int *degmixmat) {
  Vertex f, m;
  Edge e;
  Vertex *id, *od;
  Edge i, j, ne = *nedge;
  int time;
  int bipartite = *nfem;
  double *heads, *tails;
  TreeNode *ie, *oe;
  Network nw;

  /* R's serialization of matrixes is column-major, so this works: */
  nw = NetworkInitialize(edge, edge+*nedge, ne, *nnodes, 0, bipartite, 1);
  ie=nw.inedges;
  oe=nw.outedges;
  id=nw.indegree;
  od=nw.outdegree;
  
  /* Step through time one click at a time */
  for (time=j=0; time <= *ntimestep; time++) {
    /* Update the DEGMIXMAT */
    for (i=0; i < *nfem; i++) {
      switch (od[i+1]) {
        case 0: ++DEGMIXMAT(i, 0);
        break;
        case 1: /*if (id[oe[i+1].value] == 1) { */
          m = oe[i+1].value;
          e=EdgetreeMinimum(oe, i+1);
          if (id[m]==1) {
            ++DEGMIXMAT(i, 1);
          } else {
            ++DEGMIXMAT(i, 2);
          }
          break;
          default: ++DEGMIXMAT(i, 3);
      }
    }
    for (; i < *nnodes; i++) {
      switch (id[i+1]) {
        case 0: ++DEGMIXMAT(i, 0);
        break;
        case 1: if (od[ie[i+1].value] == 1) {
          ++DEGMIXMAT(i, 1);
        } else {
          ++DEGMIXMAT(i, 2);
        }
        break;
        default: ++DEGMIXMAT(i, 3);
      }
    }
    /* Toggle the edges at this timestep */
    if (time < *ntimestep) {
      for(; CHANGE(j,0) == time; j++) {
        ToggleEdgeWithTimestamp(CHANGE(j, 1), CHANGE(j, 2), &nw); 
      }
    }
  }
  /* Free memory used by network object before returning */  
  NetworkDestroy (&nw);
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
void godfather_wrapper (int *heads, int *tails, int *dnedges,
			int *dn, int *dflag, int *bipartite, 
			int *nterms, char **funnames,
			char **sonames, 
			int *totalntoggles, int *timestamps, 
			int *toggleheads, int *toggletails,
			double *inputs, 
			double *changestats, 
			int *newnetworkheads, 
			int *newnetworktails, 
			int *accumulate, 
			int *fVerbose, 
			int *maxedges) {
  int directed_flag, maxtoggles, ntoggles;
  Vertex n_nodes, bip;
  Edge i, j, n_edges, nmax, samplesize, nextedge, tnt;
  Network nw;
  Model *m;
  int thistime; 
  
  n_nodes = (Vertex)*dn; /* coerce double *dn to type Vertex */
  n_edges = (Edge)*dnedges; /* coerce double *dnedges to type Vertex */
  nmax = (Edge)*maxedges; /* coerce double *maxedges to type Vertex */
  bip = (Vertex)*bipartite; /* coerce double *bipartite to type Vertex */    
  tnt = (Edge) *totalntoggles; /* coerce double *totalntoggles to type Edge */ 
  directed_flag = *dflag;

  m=ModelInitialize(*funnames, *sonames, inputs, *nterms);
  nw = NetworkInitialize(heads, tails, n_edges, n_nodes, directed_flag, bip, 1);

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
  
  /* Now start obtaining change statistics */
  thistime = timestamps[0];
  for(i = ntoggles = 0; i <= tnt; i++){
    if (i<tnt && timestamps[i] == thistime) {
      ++ntoggles; /* ntoggles counts how many timestamps match the current value. */
    } else { /* we have encountered a new timestamp value. */
      
      /* Calculate the change statistics relative to previous step */
      changestats += m->n_stats;
      ChangeStats(ntoggles, toggleheads+i-ntoggles, toggletails+i-ntoggles, &nw, m);

      /* Accumulate change statistics */
      for (j=0; j<m->n_stats; j++){
        changestats[j] = changestats[j-m->n_stats] + m->workspace[j];
      }

      /* Make proposed toggles (for real this time) */
      for (j=ntoggles; j > 0; j--){
	if (!(*accumulate) || EdgetreeSearch(toggleheads[i-j], toggletails[i-j], nw.outedges) == 0) { 
	  ToggleEdgeWithTimestamp(toggleheads[i-j], toggletails[i-j], &nw);
	}
      }
      
      /* Finished with this timestamp; go on to next one */
      if (i<tnt) {
        ntoggles = 1;
        thistime = timestamps[i];
      }
    }
  }
  if (nmax>0) {
    /* record new generated network to pass back to R */
    newnetworkheads[0]=newnetworktails[0]=EdgeTree2EdgeList(newnetworkheads+1,newnetworktails+1,&nw,nmax);
    if (newnetworkheads[0]>=nmax) { 
      Rprintf("Error!  The value of maxedges was not set high enough in ergm.godfather\n");
    }
  }
  /* Clean up and return */
  ModelDestroy(m);
  NetworkDestroy(&nw);
}
