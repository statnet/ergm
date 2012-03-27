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
   The value maxNumDyadTypes is the largest allowable number of
   unique sets of change statistics.
 Re-rewritten by Pavel Krivitsky to make compression fast. ;)
 *****************/

/* *** don't forget tail -> head, and so this function accepts
   tails before heads now */

void MPLE_wrapper (int *tails, int *heads, int *dnedges,
		   int *dn, int *dflag, int *bipartite, int *nterms, 
		   char **funnames, char **sonames, double *inputs,  
		   int *responsevec, double *covmat,
		   int *weightsvector,
		   double * offset, double * compressedOffset,
		   int *maxNumDyadTypes) {
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
  
  MpleInit_hash(responsevec, covmat, weightsvector, offset, 
		compressedOffset, *maxNumDyadTypes, nw, m); 
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
				    int response, double offset){
  /* Cast all pointers to unsigned char pointers, since data need to 
     be fed to the hash function one byte at a time. */
  unsigned char *cnewRow = (unsigned char *) newRow,
    *cresponse = (unsigned char *) &response,
    *coffset = (unsigned char *) &offset;
  unsigned int crowLength = rowLength * sizeof(double);
  
  unsigned int hash=0;

#define HASH_LOOP(hash, keybyte){ hash+=keybyte; hash+=(hash<<10); hash^= (hash>>6); }
  for(unsigned int i=0; i<crowLength; i++) HASH_LOOP(hash, cnewRow[i]);
  for(unsigned int i=0; i<sizeof(int); i++) HASH_LOOP(hash, cresponse[i]);
  for(unsigned int i=0; i<sizeof(double); i++) HASH_LOOP(hash, coffset[i]);
#undef HASH_LOOP

  hash += (hash<<3);
  hash ^= (hash>>11);
  hash += (hash<<15);

  return(hash % numRows);
}

/*R_INLINE*/ unsigned int insCovMatRow(double *newRow, double *matrix, unsigned int rowLength, unsigned int numRows,
			  int response, int *responsevec,
			  double offset, double *compressedOffset, int *weights ){
  unsigned int hash_pos = hashCovMatRow(newRow, rowLength, numRows, response, offset), pos, round;
  
  for(/*unsigned int*/ pos=hash_pos, round=0; !round ; pos = (pos+1)%numRows, round+=(pos==hash_pos)?1:0){
//    Rprintf("pos %d round %d hash_pos %d\n",pos,round,hash_pos);
    if(weights[pos]==0){ /* Space is unoccupied. */
      weights[pos]=1;
      compressedOffset[pos]=offset;
      responsevec[pos]=response;
      memcpy(matrix+rowLength*pos,newRow,rowLength*sizeof(double));
      return TRUE;
    }else {
      
      if( compressedOffset[pos]==offset &&
	      responsevec[pos]==response &&
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
*****************/


void MpleInit_hash(int *responsevec, double *covmat, int *weightsvector,
		   double *offset, double *compressedOffset,
		   int maxNumDyadTypes, Network *nwp, Model *m) {
  Edge dyadNum=0;
  Vertex rowmax;
  ModelTerm *mtp;
  double *newRow = (double *) R_alloc(m->n_stats,sizeof(double));
  /* Note:  This function uses macros found in changestat.h */
  
  for(Vertex t = 1; t <= (BIPARTITE ? BIPARTITE : N_NODES); t++){
    for(Vertex h = MAX(DIRECTED ? 0 : t, BIPARTITE)+1; h <= N_NODES; h++){
      if(h==t) continue;
      int response = IS_OUTEDGE(t,h);
      unsigned int totalStats = 0;
      /* Let mtp loop through each model term */
      for (mtp=m->termarray; mtp < m->termarray + m->n_terms; mtp++){
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
		       maxNumDyadTypes, response, 
		       responsevec, offset ? offset[dyadNum++]:0, 
		       compressedOffset, weightsvector)) {
	error("Too many unique dyads!");
      }
    }
  }
}
