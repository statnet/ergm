/*  File src/allstatistics.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2017 Statnet Commons
 */
#include "MPLE.h"
#include "changestat.h"

void RecurseOffOn(int *nodelist1,int *nodelist2, int nodelistlength, 
       int currentnodes, double *changeStats, double *cumulativeStats,
       double *covmat, int *weightsvector,
       int maxNumDyadTypes, Network *nwp, Model *m);

unsigned int InsNetStatRow(double *newRow, double *matrix, int rowLength, 
       int numRows, int *weights );

unsigned int hashNetStatRow(double *newRow, int rowLength, 
       int numRows);

        
/* *****************
 void AllStatistics

 Wrapper for a call from R.  Based on MPLE_wrapper but much simpler.  
 Produces matrix of network statistics for an arbitrary statnet model
 by an algorithm that starts with the network passed in the
 formula, then recursively toggling each edge two times so that every 
 possible network is visited.
 *****************/

void AllStatistics (
       int *tails, 
       int *heads,
       int *dnedges,
		   int *dn, /* Number of nodes */
       int *dflag, /* directed flag */
       int *bipartite,
       int *nterms, 
		   char **funnames, 
       char **sonames, 
       double *inputs,  
       double *covmat,
		   int *weightsvector,
       int *maxNumDyadTypes) {
  Network nw, *nwp;

  Vertex n_nodes = (Vertex) *dn; 
  int directed_flag = *dflag;
  int nodelistlength, rowmax, *nodelist1, *nodelist2;
  Vertex bip = (Vertex) *bipartite;
  Model *m;
  ModelTerm *mtp;

  /* Step 1:  Initialize empty network and initialize model */
  GetRNGstate(); /* Necessary for R random number generator */
  nw=NetworkInitialize(tails, heads, *dnedges,
		       n_nodes, directed_flag, bip, 0, 0, NULL);
  nwp = &nw;
  m=ModelInitialize(*funnames, *sonames, &inputs, *nterms);
  
  /* Step 2:  Build nodelist1 and nodelist2, which together give all of the
  dyads in the network. */
  if (BIPARTITE > 0) { /* Assuming undirected in the bipartite case */
    nodelistlength = BIPARTITE * (N_NODES-BIPARTITE);
    rowmax = BIPARTITE + 1;
  } else {
    nodelistlength = N_NODES * (N_NODES-1) / (DIRECTED? 1 : 2);
    rowmax = N_NODES;
  }
  nodelist1 = (int *) R_alloc(nodelistlength, sizeof(int));
  nodelist2 = (int *) R_alloc(nodelistlength, sizeof(int));
  int count = 0;
  for(int i=1; i < rowmax; i++) {
    for(int j = MAX(i,BIPARTITE)+1; j <= N_NODES; j++) {
      for(int d=0; d <= DIRECTED; d++) { /*trivial loop if undirected*/
        nodelist1[count] = d==1? j : i;
        nodelist2[count] = d==1? i : j;
        count++;
      }
    }
  }

  /* Step 3:  Initialize values of mtp->dstats so they point to the correct
  spots in the newRow vector.  These values will never change. */
  double *changeStats     = (double *) R_alloc(m->n_stats,sizeof(double));
  double *cumulativeStats = (double *) R_alloc(m->n_stats,sizeof(double));
  for (int i=0; i < m->n_stats; i++) cumulativeStats[i]=0.0;

  unsigned int totalStats = 0;
  for (mtp=m->termarray; mtp < m->termarray + m->n_terms; mtp++){
    mtp->dstats = changeStats + totalStats;
    /* Update mtp->dstats pointer to skip atail by mtp->nstats */
    totalStats += mtp->nstats; 
  }
  if (totalStats != m->n_stats) {
    Rprintf("I thought totalStats=%d and m->nstats=%d should be the same.\n", 
    totalStats, m->n_stats);
  }

  /* Step 4:  Begin recursion */
  RecurseOffOn(nodelist1, nodelist2, nodelistlength, 0, changeStats, 
           cumulativeStats, covmat, weightsvector, *maxNumDyadTypes, nwp, m);

  /* Step 5:  Deallocate memory and return */
  ModelDestroy(m);
  NetworkDestroy(nwp);
  PutRNGstate(); /* Must be called after GetRNGstate before returning to R */
}

void RecurseOffOn(
       int *nodelist1, 
       int *nodelist2, 
       int nodelistlength, 
       int currentnodes, 
       double *changeStats,
       double *cumulativeStats,
       double *covmat,
       int *weightsvector,
		   int maxNumDyadTypes, 
       Network *nwp, 
       Model *m) {
  ModelTerm *mtp;

  /* Loop through twice for each dyad: Once for edge and once for no edge */
  for (int i=0; i<2; i++) {
    /* Recurse if currentnodes+1 is not yet nodelistlength */
    if (currentnodes+1 < nodelistlength) {
      RecurseOffOn(nodelist1, nodelist2, nodelistlength, currentnodes+1, changeStats, 
              cumulativeStats, covmat, weightsvector, maxNumDyadTypes, nwp, m);
    } else { /* Add newRow of statistics to hash table */
      //Rprintf("Here's network number %d with changestat %f\n",counter++, changeStats[0]);
      //NetworkEdgeList(nwp);
      if(!InsNetStatRow(cumulativeStats, covmat, m->n_stats, maxNumDyadTypes, weightsvector)) {
        error("Too many unique dyads!");
      }
    }

    /* Calculate the change statistic(s) associated with toggling the 
       dyad represented by nodelist1[currentnodes], nodelist2[currentnodes] */
    for (mtp=m->termarray; mtp < m->termarray + m->n_terms; mtp++){
      (*(mtp->d_func))(1, nodelist1+currentnodes, nodelist2+currentnodes, mtp, nwp);
    }
    for (int j=0; j < m->n_stats; j++) cumulativeStats[j] += changeStats[j];
    /* Now toggle the dyad so it's ready for the next pass */
    ToggleEdge(nodelist1[currentnodes], nodelist2[currentnodes], nwp);
  }
}

unsigned int InsNetStatRow(
                     double *newRow, 
                     double *matrix, 
                     int rowLength, 
                     int numRows, 
                     int *weights ){
  unsigned int pos, round;
  unsigned int hash_pos = hashNetStatRow(newRow, rowLength, numRows);
  
  for(pos=hash_pos, round=0; !round ; pos = (pos+1)%numRows, round+=(pos==hash_pos)?1:0){
    if(weights[pos]==0){ /* Space is unoccupied. */
      weights[pos]=1;
      memcpy(matrix+rowLength*pos,newRow,rowLength*sizeof(double));
      return TRUE;
    }else{
      if(memcmp(matrix+rowLength*pos, newRow, rowLength*sizeof(double))==0 ){ /* Rows are identical. */
        weights[pos]++;
        return TRUE;
      }
    }
  }
  return FALSE; /* Insertion unsuccessful: the table is full. */
}


/*************
Hashes the covariates, offset, and response onto an unsigned integer in the interval [0,numRows).
Uses Jenkins One-at-a-Time hash.

numRows should, ideally, be a power of 2, but it doesn't have to be.
**************/
unsigned int hashNetStatRow(double *newRow, int rowLength, int numRows) {
  /* Cast all pointers to unsigned char pointers, since data need to 
     be fed to the hash function one byte at a time. */
  unsigned char *cnewRow = (unsigned char *) newRow;
  unsigned int crowLength = rowLength * sizeof(double);
  unsigned int hash=0;

#define HASH_LOOP(hash, keybyte){ hash+=keybyte; hash+=(hash<<10); hash^= (hash>>6); }
  for(unsigned int i=0; i<crowLength; i++) HASH_LOOP(hash, cnewRow[i]);
#undef HASH_LOOP

  hash += (hash<<3);
  hash ^= (hash>>11);
  hash += (hash<<15);

  return(hash % numRows);
}



