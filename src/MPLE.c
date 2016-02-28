/*  File src/MPLE.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2013 Statnet Commons
 */
#include "MPLE.h"
#include "changestat.h"

/* *****************
 void MPLE_wrapper

 Wrapper for a call from R.  
 Only find the MPLE, so no MCMC stuff is necessary. 

 Re-Written by David Schruth to add compression:
   Since there are often many different 0->1 dyad changes that
   result in identical change statistic vectors, this routine
   now returns only the unique change statistic vectors in the 
   covmat matrix, along with a separate vector, weightsvector,
   giving the number of repetitions of each unique set of change
   statistics.  Note that two sets of statistics are considered
   unique in this context if the corresponding response values
   (i.e., dyad values, edge or no edge) are unequal.
   The value maxDyadTypes is the largest allowable number of
   unique sets of change statistics.
 Re-rewritten by Pavel Krivitsky to make compression fast. ;)
 *****************/

/* *** don't forget tail -> head, and so this function accepts
   tails before heads now */

void MPLE_wrapper(int *tails, int *heads, int *dnedges,
		  int *wl, int *lel,
		  int *dn, int *dflag, int *bipartite, int *nterms, 
		  char **funnames, char **sonames, double *inputs,  
		  int *responsevec, double *covmat,
		  int *weightsvector,
		  int *maxDyads, int *maxDyadTypes){
  Network nw[2];
  Vertex n_nodes = (Vertex) *dn; 
  Edge n_edges = (Edge) *dnedges;
  int directed_flag = *dflag;
  Vertex bip = (Vertex) *bipartite;
  Model *m;

  GetRNGstate(); /* Necessary for R random number generator */
  nw[0]=NetworkInitialize(tails, heads, n_edges,
                          n_nodes, directed_flag, bip, 0, 0, NULL);
  m=ModelInitialize(*funnames, *sonames, &inputs, *nterms);
  
  if(*wl) MpleInit_hash_wl(responsevec, covmat, weightsvector, lel, *maxDyadTypes, nw, m); 
  else MpleInit_hash_bl(responsevec, covmat, weightsvector, lel, *maxDyads, *maxDyadTypes, nw, m); 

  ModelDestroy(m);
  NetworkDestroy(nw);
  PutRNGstate(); /* Must be called after GetRNGstate before returning to R */
}

/*************
Hashes the covariates, offset, and response onto an unsigned integer in the interval [0,numRows).
Uses Jenkins One-at-a-Time hash.

numRows should, ideally, be a power of 2, but doesn't have to be.
**************/
/*R_INLINE*/ unsigned int hashCovMatRow(double *newRow, unsigned int rowLength, unsigned int numRows,
				    int response){
  /* Cast all pointers to unsigned char pointers, since data need to 
     be fed to the hash function one byte at a time. */
  unsigned char *cnewRow = (unsigned char *) newRow,
    *cresponse = (unsigned char *) &response;
  unsigned int crowLength = rowLength * sizeof(double);
  
  unsigned int hash=0;

#define HASH_LOOP(hash, keybyte){ hash+=keybyte; hash+=(hash<<10); hash^= (hash>>6); }
  for(unsigned int i=0; i<crowLength; i++) HASH_LOOP(hash, cnewRow[i]);
  for(unsigned int i=0; i<sizeof(int); i++) HASH_LOOP(hash, cresponse[i]);
#undef HASH_LOOP

  hash += (hash<<3);
  hash ^= (hash>>11);
  hash += (hash<<15);

  return(hash % numRows);
}

/*R_INLINE*/ unsigned int insCovMatRow(double *newRow, double *matrix, unsigned int rowLength, unsigned int numRows,
			  int response, int *responsevec,
			  int *weights ){
  unsigned int hash_pos = hashCovMatRow(newRow, rowLength, numRows, response), pos, round;
  
  for(/*unsigned int*/ pos=hash_pos, round=0; !round ; pos = (pos+1)%numRows, round+=(pos==hash_pos)?1:0){
//    Rprintf("pos %d round %d hash_pos %d\n",pos,round,hash_pos);
    if(weights[pos]==0){ /* Space is unoccupied. */
      weights[pos]=1;
      responsevec[pos]=response;
      memcpy(matrix+rowLength*pos,newRow,rowLength*sizeof(double));
      return TRUE;
    }else {
      
      if(responsevec[pos]==response &&
      memcmp(matrix+rowLength*pos,newRow,rowLength*sizeof(double))==0 ){ /* Rows are identical. */
        weights[pos]++;
        return TRUE;
      }
    }
  }
  return FALSE; /* Insertion unsuccessful: the table is full. */
}

/*****************
 void MpleInit_hash

 For finding the MPLE, an extra bit of initialization is required:  
 we must build the matrix of covariates to be used in the logistic 
 regression routine.  This matrix has #rows equal to the number of 
 possible edges in the network (n choose 2 for an undirected network) 
 and #cols equal to the number of network statistics in the model, 
 which is also the number of parameters.  The row for each edge in 
 the matrix should contain the changes in each of the statistics that 
 would occur if the given edge is toggled from 0 to 1, leaving all 
 other edges as they are in the observed network.  The response vector 
 for the logistic regression is simply the vector of indicators 
 giving the states of the edges in the observed network.

 It comes in two versions: one that iterates through all valid dyads,
 skipping those on the additional list, which is treated as a
 blacklist; and another that iterates only through dyads on the
 additional ist, which is treated as a whitelist.

*****************/

/* Maps a dyad ID to a tail and a head, serializing them in a
   tail-major order.  Here, __s represent structural 0s.

   Note that both undirected and directed unipartite networks are
   serialized in the exact same manner. The t<h "filtering" or
   swapping must be done later. (The alternative appears to require a
   square root, which is inefficient.)

   Bipartite(2 b1s, 3 b2s):

          h

   __ __ 00 01 02
t  __ __ 03 04 05
   __ __ __ __ __
   __ __ __ __ __
   __ __ __ __ __

   Unipartite:

          h

   __ 00 01 02 03
   04 __ 05 06 07
t  08 09 __ 10 11 
   12 13 14 __ 15 
   16 17 18 19 __    

 */

// Note that id MUST be an integer of some sort, or you'll get floating point division.
#define DYADID2T(id) (BIPARTITE ? id / (N_NODES-BIPARTITE) + 1 : id / (N_NODES-1) + 1)
#define DYADID2H(id) (BIPARTITE ? id % (N_NODES-BIPARTITE) + 1 + BIPARTITE : \
		      id % (N_NODES-1) + 1 + (id % (N_NODES-1) >= id / (N_NODES-1)))


// Euclid's Algorithm to compute Greatest Common Divisor of a and b.
unsigned long int gcd(unsigned long int a, unsigned long int b){
  if(b==0) return a;
  else return gcd(b, a%b);
}

void MpleInit_hash_bl(int *responsevec, double *covmat, int *weightsvector,
		      int *bl, 
		      Edge maxDyads, Edge maxDyadTypes, Network *nwp, Model *m){
  double *newRow = (double *) R_alloc(m->n_stats,sizeof(double));
  /* Note:  This function uses macros found in changestat.h */

  Dyad dc = DYADCOUNT(N_NODES, BIPARTITE, DIRECTED);
  if(!BIPARTITE && !DIRECTED){
    dc*=2;
    maxDyads*=2;
  }
  
  // Find a number relatively prime with the dyad count:
  Edge step = MAX(N_NODES/3,2);
  while(gcd(dc,step)!=1) step++;

  for(Dyad d = 0, di = 0; d < MIN(maxDyads,dc); d++, di=(di+step)%dc){
    Vertex t=DYADID2T(di), h=DYADID2H(di);
    //Rprintf("d=%d di=%d t=%d h=%d ", d,di,t,h);
    if((!DIRECTED && t>=h) || iEdgeListSearch(t,h,bl)) continue;

    int response = IS_OUTEDGE(t,h);
    unsigned int totalStats = 0;
    /* Let mtp loop through each model term */
    for (ModelTerm *mtp=m->termarray; mtp < m->termarray + m->n_terms; mtp++){
      mtp->dstats = newRow + totalStats;
      /* Now call d_xxx function, which updates mtp->dstats to reflect
	 changing the current dyad.  */
      (*(mtp->d_func))(1, &t, &h, mtp, nwp);
      /* dstats values reflect changes in current dyad; for MPLE, 
	 values must reflect going from 0 to 1.  Thus, we have to reverse 
	 the sign of dstats whenever the current edge exists. */
      if(response){
	  for(unsigned int l=0; l<mtp->nstats; l++){
	    mtp->dstats[l] = -mtp->dstats[l];
	  }
      }
      /* Update mtp->dstats pointer to skip ahead by mtp->nstats */
      totalStats += mtp->nstats; 
    }
    /* In next command, if there is an offset vector then its total
       number of entries should match the number of times through the 
       inner loop (i.e., the number of dyads in the network) */          
    if(!insCovMatRow(newRow, covmat, m->n_stats,
		     maxDyadTypes, response, 
		     responsevec, weightsvector)) {
      warning("Too many unique dyads. MPLE is approximate, and MPLE standard errors are suspect.");
      break;
    }
  }
}

void MpleInit_hash_wl(int *responsevec, double *covmat, int *weightsvector,
		      int *wl, 
		      Edge maxDyadTypes, Network *nwp, Model *m){
  double *newRow = (double *) R_alloc(m->n_stats,sizeof(double));
  /* Note:  This function uses macros found in changestat.h */
  Edge wledges = wl[0];
  wl++;

  for(Edge wlpos = 0; wlpos < wledges; wlpos++){
    Vertex t=wl[wlpos], h=(wl+wledges)[wlpos];
    
    int response = IS_OUTEDGE(t,h);
    unsigned int totalStats = 0;
    /* Let mtp loop through each model term */
    for (ModelTerm *mtp=m->termarray; mtp < m->termarray + m->n_terms; mtp++){
      mtp->dstats = newRow + totalStats;
      /* Now call d_xxx function, which updates mtp->dstats to reflect
	 changing the current dyad.  */
      (*(mtp->d_func))(1, &t, &h, mtp, nwp);
      /* dstats values reflect changes in current dyad; for MPLE, 
	 values must reflect going from 0 to 1.  Thus, we have to reverse 
	 the sign of dstats whenever the current edge exists. */
      if(response){
	for(unsigned int l=0; l<mtp->nstats; l++){
	  mtp->dstats[l] = -mtp->dstats[l];
	}
      }
      /* Update mtp->dstats pointer to skip ahead by mtp->nstats */
      totalStats += mtp->nstats; 
    }
    /* In next command, if there is an offset vector then its total
       number of entries should match the number of times through the 
       inner loop (i.e., the number of dyads in the network) */          
    if(!insCovMatRow(newRow, covmat, m->n_stats,
		     maxDyadTypes, response, 
		     responsevec, weightsvector)) {
      warning("Too many unique dyads. MPLE is approximate, and MPLE standard errors are suspect.");
      break;
    }
  }
}
