/*
 *  File ergm/src/MPLE.c
 *  Part of the statnet package, http://statnetproject.org
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) in
 *    http://statnetproject.org/attribution
 *
 * Copyright 2003 Mark S. Handcock, University of Washington
 *                David R. Hunter, Penn State University
 *                Carter T. Butts, University of California - Irvine
 *                Steven M. Goodreau, University of Washington
 *                Martina Morris, University of Washington
 * Copyright 2007 The statnet Development Team
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
   The value maxNumDyadTypes is the largest allowable number of
   unique sets of change statistics.
 Re-rewritten by Pavel Krivitsky to make compression fast. ;)
 *****************/

void MPLE_wrapper (int *heads, int *tails, int *dnedges,
       int *maxpossibleedges,
		   int *dn, int *dflag, int *bipartite, int *nterms, 
		   char **funnames, char **sonames, double *inputs,  
		   int *responsevec, double *covmat,
		   int *weightsvector,
		   double * offset, double * compressedOffset,
		   int *maxNumDyadTypes, int *maxMPLEsamplesize, int *compressflag) {
  Network nw[2];
  Vertex n_nodes = (Vertex) *dn; 
  Edge n_edges = (Edge) *dnedges;
  int directed_flag = *dflag;
  int hammingterm;
  Vertex bip = (Vertex) *bipartite;
  Edge maxMPLE = (Edge) *maxMPLEsamplesize;
  Vertex hhead, htail;
  Edge  nddyads, kedge;
  Model *m;
  ModelTerm *thisterm;

  GetRNGstate(); /* Necessary for R random number generator */
  nw[0]=NetworkInitialize(heads, tails, n_edges,
                          n_nodes, directed_flag, bip, 0);
  m=ModelInitialize(*funnames, *sonames, inputs, *nterms);
  
  hammingterm=ModelTermHamming (*funnames, *nterms);
  if(hammingterm>0){
   Network nwhamming;
   thisterm = m->termarray + hammingterm - 1;
   nddyads = (Edge)(thisterm->inputparams[0]);
   nwhamming=NetworkInitializeD(thisterm->inputparams+1, 
				thisterm->inputparams+1+nddyads,
			       	nddyads, n_nodes, directed_flag, bip,0);
   nddyads=0;
   nw[1]=NetworkInitializeD(thisterm->inputparams+1, 
			   thisterm->inputparams+1+nddyads, nddyads,
         n_nodes, directed_flag, bip,0);
   for (kedge=1; kedge <= nwhamming.nedges; kedge++) {
     FindithEdge(&hhead, &htail, kedge, &nwhamming);
     if(EdgetreeSearch(hhead, htail, nw[0].outedges) == 0){
       ToggleEdge(hhead, htail, &nw[1]);
     }
   }
   for (kedge=1; kedge <= nw[0].nedges; kedge++) {
     FindithEdge(&hhead, &htail, kedge, &nw[0]);
     if(EdgetreeSearch(hhead, htail, nwhamming.outedges) == 0){
       ToggleEdge(hhead, htail, &nw[1]);
     }
   }
/*   Rprintf("Initial number of discordant %d Number of g0 ties %d Number of ties in g %d\n",nw[1].nedges, nwhamming.nedges,nw[0].nedges); */
   NetworkDestroy(&nwhamming);
  }

  if (*compressflag) 
    MpleInit_hash(responsevec, covmat, weightsvector, offset, 
		  compressedOffset, *maxNumDyadTypes, maxMPLE, nw, m); 
  else
    MpleInit_no_compress(responsevec, covmat, weightsvector, 
      *maxNumDyadTypes, maxMPLE, nw, m);
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
 void MpleInit_*

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

 The *_hash version also "compresses" the output by tabulating
 duplicate rows.
*****************/

void MpleInit_no_compress (int *responsevec, double *covmat, int *weightsvector,
		     int maxNumDyadTypes, Edge maxMPLE, Network *nwp, Model *m) {
  int l, d, outflag = 0, inflag = 0, thisRowNumber, totalStats, *currentResponse;
  double *covMatPosition;
  Vertex i, j , rowmax;
  ModelTerm *mtp;
  /* Note:  This function uses macros found in changestats.h */
  
  covMatPosition = covmat;
  currentResponse = responsevec;
  thisRowNumber = 0;
  if(BIPARTITE > 0) rowmax = BIPARTITE + 1;
  else              rowmax = N_NODES;
  for(i=1; i < rowmax; i++){
    for(j = MAX(i,BIPARTITE)+1; j <= N_NODES; j++){
      for(d=0; d <= DIRECTED; d++){ /*trivial loop if undirected*/
        if (d==1)*currentResponse = inflag = IS_INEDGE(i,j);
        else     *currentResponse = outflag = IS_OUTEDGE(i,j);
        totalStats = 0;
        if(*currentResponse || i <= maxMPLE){   
          /* Let mtp loop through each model term */
          for (mtp=m->termarray; mtp < m->termarray + m->n_terms; mtp++){
            mtp->dstats = covMatPosition + totalStats;
            /* Now call d_xxx function, which updates mtp->dstats to reflect
            changing the current dyad.  */
            if(d==0) (*(mtp->d_func))(1, &i, &j, mtp, nwp);
            else(*(mtp->d_func))(1, &j, &i, mtp, nwp);
            /* dstats values reflect changes in current dyad; for MPLE, 
            values must reflect going from 0 to 1.  Thus, we have to reverse 
            the sign of dstats whenever the current edge exists. */
            if((d==0 && outflag) || (d==1 && inflag)){
              for(l=0; l<mtp->nstats; l++){
                mtp->dstats[l] = -mtp->dstats[l];
              }
            }
            /* Update mtp->dstats pointer to skip ahead by mtp->nstats */
            totalStats += mtp->nstats; 
          }
          if(thisRowNumber<maxNumDyadTypes){ 
            /* Shift the pointer n parameters forward in
            the covariate matrix vector */
            covMatPosition += m->n_stats; /* New row in covmat matrix */
            currentResponse++; /* New response value */
            weightsvector[thisRowNumber++]=1; /* New # unique rows */
          } else{ /* Do nothing for now if thisRowNumber >=maxNumDyadTypes */ 
          }
        }
      }
    }
  }
}

void MpleInit_hash(int *responsevec, double *covmat, int *weightsvector,
		   double *offset, double *compressedOffset,
		   int maxNumDyadTypes, Edge maxMPLE, Network *nwp, Model *m) {
  int outflag = 0, inflag = 0;
  Edge dyadNum=0;
  Vertex rowmax;
  ModelTerm *mtp;
  double *newRow = (double *) R_alloc(m->n_stats,sizeof(double));
  /* Note:  This function uses macros found in changestats.h */
  
  if(BIPARTITE > 0) rowmax = BIPARTITE + 1;
  else              rowmax = N_NODES;
  for(Vertex i=1; i < rowmax; i++){
    for(Vertex j = MAX(i,BIPARTITE)+1; j <= N_NODES; j++){
      for(unsigned int d=0; d <= DIRECTED; d++){ /*trivial loop if undirected*/
        int response;
        if (d==1) response = inflag = IS_INEDGE(i,j);
        else      response = outflag = IS_OUTEDGE(i,j);
        unsigned int totalStats = 0;
        if(response || i <= maxMPLE){   
          /* Let mtp loop through each model term */
          for (mtp=m->termarray; mtp < m->termarray + m->n_terms; mtp++){
            mtp->dstats = newRow + totalStats;
            /* Now call d_xxx function, which updates mtp->dstats to reflect
            changing the current dyad.  */
            if(d==0) (*(mtp->d_func))(1, &i, &j, mtp, nwp);
            else(*(mtp->d_func))(1, &j, &i, mtp, nwp);
            /* dstats values reflect changes in current dyad; for MPLE, 
            values must reflect going from 0 to 1.  Thus, we have to reverse 
            the sign of dstats whenever the current edge exists. */
            if((d==0 && outflag) || (d==1 && inflag)){
              for(unsigned int l=0; l<mtp->nstats; l++){
                mtp->dstats[l] = -mtp->dstats[l];
              }
            }
            /* Update mtp->dstats pointer to skip ahead by mtp->nstats */
            totalStats += mtp->nstats; 
          }
          if(!insCovMatRow(newRow, covmat, m->n_stats,
			   maxNumDyadTypes, response, 
			   responsevec, offset ? offset[dyadNum++]:0, 
			   compressedOffset, weightsvector)) {
            error("Too many unique dyads!");
          }
        }
      }
    }
  }
}
