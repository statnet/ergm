/*  File src/allstatistics.c in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2021 Statnet Commons
 */
#include "MPLE.h"
#include "ergm_changestat.h"
#include "ergm_state.h"
#include "ergm_util.h"

void RecurseOffOn(ErgmState *s, Vertex *nodelist1,Vertex *nodelist2, Vertex nodelistlength, 
       Vertex currentnodes, double *changeStats, double *cumulativeStats,
       double *covmat, unsigned int *weightsvector,
       unsigned int maxNumDyadTypes);

unsigned int InsNetStatRow(double *newRow, double *matrix, unsigned int rowLength, 
       unsigned int numRows, unsigned int *weights );

unsigned int hashNetStatRow(double *newRow, unsigned int rowLength, 
       unsigned int numRows);

        
/* *****************
 void AllStatistics

 Wrapper for a call from R.  Based on MPLE_wrapper but much simpler.  
 Produces matrix of network statistics for an arbitrary statnet model
 by an algorithm that starts with the network passed in the
 formula, then recursively toggling each edge two times so that every 
 possible network is visited.
 *****************/

SEXP AllStatistics(SEXP stateR,
                   // Allstats settings
                   SEXP maxNumDyadTypes){

  /* Step 1:  Initialize empty network and initialize model */
  GetRNGstate(); /* Necessary for R random number generator */
  ErgmState *s = ErgmStateInit(stateR, ERGM_STATE_NO_INIT_PROP);

  Network *nwp = s->nwp;
  Model *m = s->m;

  /* Step 2:  Build nodelist1 and nodelist2, which together give all of the
  dyads in the network. */
  Dyad nodelistlength = DYADCOUNT(nwp);
  Vertex rowmax;
  if (BIPARTITE > 0) { /* Assuming undirected in the bipartite case */
    rowmax = BIPARTITE + 1;
  } else {
    rowmax = N_NODES;
  }
  Vertex *nodelist1 = (Vertex *) R_alloc(nodelistlength, sizeof(int));
  Vertex *nodelist2 = (Vertex *) R_alloc(nodelistlength, sizeof(int));
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
  EXEC_THROUGH_TERMS(m, {
    mtp->dstats = changeStats + totalStats;
    /* Update mtp->dstats pointer to skip atail by mtp->nstats */
    totalStats += mtp->nstats; 
    });
  if (totalStats != m->n_stats) {
    Rprintf("I thought totalStats=%d and m->n_stats=%d should be the same.\n", 
    totalStats, m->n_stats);
  }

  SEXP stats = PROTECT(allocVector(REALSXP, asInteger(maxNumDyadTypes)*m->n_stats));
  memset(REAL(stats), 0, asInteger(maxNumDyadTypes)*m->n_stats*sizeof(double));
  SEXP weights = PROTECT(allocVector(INTSXP, asInteger(maxNumDyadTypes)));
  memset(INTEGER(weights), 0, asInteger(maxNumDyadTypes)*sizeof(int));
  SEXP outl = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(outl, 0, stats);
  SET_VECTOR_ELT(outl, 1, weights);

  /* Step 4:  Begin recursion */
  RecurseOffOn(s, nodelist1, nodelist2, nodelistlength, 0, changeStats, 
	       cumulativeStats, REAL(stats), (unsigned int*) INTEGER(weights), asInteger(maxNumDyadTypes));

  /* Step 5:  Deallocate memory and return */
  ErgmStateDestroy(s);
  PutRNGstate(); /* Must be called after GetRNGstate before returning to R */
  UNPROTECT(3);
  return outl;
}

static unsigned int interrupt_steps = 0; // It's OK if this overflows.

void RecurseOffOn(ErgmState *s,
       Vertex *nodelist1, 
       Vertex *nodelist2, 
       Vertex nodelistlength, 
       Vertex currentnodes, 
       double *changeStats,
       double *cumulativeStats,
       double *covmat,
       unsigned int *weightsvector,
       unsigned int maxNumDyadTypes) {

  R_CheckUserInterruptEvery(1024u, interrupt_steps++);
  Network *nwp = s->nwp;
  Model *m = s->m;

  /* Loop through twice for each dyad: Once for edge and once for no edge */
  for (int i=0; i<2; i++) {
    /* Recurse if currentnodes+1 is not yet nodelistlength */
    if (currentnodes+1 < nodelistlength) {
      RecurseOffOn(s, nodelist1, nodelist2, nodelistlength, currentnodes+1, changeStats, 
              cumulativeStats, covmat, weightsvector, maxNumDyadTypes);
    } else { /* Add newRow of statistics to hash table */
      //Rprintf("Here's network number %d with changestat %f\n",counter++, changeStats[0]);
      //NetworkEdgeList(nwp);
      if(!InsNetStatRow(cumulativeStats, covmat, m->n_stats, maxNumDyadTypes, weightsvector)) {
        error("Too many unique dyads!");
      }
    }

    /* Calculate the change statistic(s) associated with toggling the 
       dyad represented by nodelist1[currentnodes], nodelist2[currentnodes] */
    Rboolean edgestate = IS_OUTEDGE((Vertex)nodelist1[currentnodes], (Vertex)nodelist2[currentnodes]);
    EXEC_THROUGH_TERMS(m, {
	if(mtp->c_func){
	  ZERO_ALL_CHANGESTATS();
	  (*(mtp->c_func))((Vertex)nodelist1[currentnodes], (Vertex)nodelist2[currentnodes], mtp, nwp, edgestate);
	}else if(mtp->d_func){
	  (*(mtp->d_func))(1, (Vertex*)nodelist1+currentnodes, (Vertex*)nodelist2+currentnodes, mtp, nwp);
	}
      });

    addonto(cumulativeStats, changeStats, m->n_stats);
    /* Now toggle the dyad so it's ready for the next pass */
    /* Inform u_* functions that the network is about to change. */
    ToggleKnownEdge(nodelist1[currentnodes], nodelist2[currentnodes], nwp, edgestate);
  }
}

unsigned int InsNetStatRow(
                     double *newRow, 
                     double *matrix, 
                     unsigned int rowLength, 
                     unsigned int numRows, 
                     unsigned int *weights ){
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
unsigned int hashNetStatRow(double *newRow, unsigned int rowLength, unsigned int numRows) {
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



