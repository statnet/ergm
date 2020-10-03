/*  File src/MPLE.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2020 Statnet Commons
 */
#include "MPLE.h"
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
		  double *wl,
		  int *dn, int *dflag, int *bipartite, int *nterms, 
		  char **funnames, char **sonames, double *inputs,  
		  int *responsevec, double *covmat,
		  int *weightsvector,
		  int *maxDyads, int *maxDyadTypes){
  Network *nwp;
  Vertex n_nodes = (Vertex) *dn; 
  Edge n_edges = (Edge) *dnedges;
  int directed_flag = *dflag;
  Vertex bip = (Vertex) *bipartite;
  Model *m;
  double *tmp = wl;
  RLEBDM1D wlm = unpack_RLEBDM1D(&tmp, n_nodes);

  GetRNGstate(); /* Necessary for R random number generator */
  nwp=NetworkInitialize((Vertex*)tails, (Vertex*)heads, n_edges,
                          n_nodes, directed_flag, bip, 0, 0, NULL);
  m=ModelInitialize(*funnames, *sonames, &inputs, *nterms);
  
  MpleInit_hash_wl_RLE(responsevec, covmat, weightsvector, &wlm, *maxDyads, *maxDyadTypes, nwp, m); 

  ModelDestroy(m);
  NetworkDestroy(nwp);
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

static inline unsigned int insCovMatRow(double *newRow, double *matrix, unsigned int rowLength, unsigned int numRows,
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

// Euclid's Algorithm to compute Greatest Common Divisor of a and b.
Dyad gcd(Dyad a, Dyad b){
  if(b==0) return a;
  else return gcd(b, a%b);
}

void MpleInit_hash_wl_RLE(int *responsevec, double *covmat, int *weightsvector,
			  RLEBDM1D *wl, 
			  Edge maxDyads, Edge maxDyadTypes, Network *nwp, Model *m){
  double *newRow = (double *) R_alloc(m->n_stats,sizeof(double));
  /* Note:  This function uses macros found in changestat.h */

  // Number of free dyads.
  Dyad dc = wl->ndyads;
  
  // Find a number relatively prime with the free dyad count:
  Edge step = MAX(N_NODES/3,2);
  while(gcd(dc,step)!=1) step++;

  { // Start a scope for loop variables.
    Vertex t, h;
    GetRandRLEBDM1D_RS(&t,&h, wl); // Find a random starting point.
    Dyad d = TH2Dyad(nwp->nnodes, t,h);
    RLERun r=0;
    
    for(Dyad i = 0; i < MIN(maxDyads,dc); i++, d=NextRLEBDM1D(d, step, wl, &r)){
      Dyad2TH(&t, &h, d, N_NODES);
      
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
  } // End scope for loop variables.
}
