/*  File src/MPLE.c in package ergm, part of the
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
#include "ergm_rlebdm.h"

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

SEXP MPLE_wrapper(SEXP stateR,
                  // MPLE settings
                  SEXP wl,
		  SEXP maxNumDyads, SEXP maxNumDyadTypes){
  GetRNGstate(); /* Necessary for R random number generator */
  ErgmState *s = ErgmStateInit(stateR, ERGM_STATE_NO_INIT_PROP);

  Model *m = s->m;

  double *tmp = REAL(wl);
  RLEBDM1D wlm = unpack_RLEBDM1D(&tmp);

  SEXP responsevec = PROTECT(allocVector(INTSXP, asInteger(maxNumDyadTypes)));
  memset(INTEGER(responsevec), 0, asInteger(maxNumDyadTypes)*sizeof(int));
  SEXP covmat = PROTECT(allocVector(REALSXP, asInteger(maxNumDyadTypes)*m->n_stats));
  memset(REAL(covmat), 0, asInteger(maxNumDyadTypes)*m->n_stats*sizeof(double));
  SEXP weightsvector = PROTECT(allocVector(INTSXP, asInteger(maxNumDyadTypes)));
  memset(INTEGER(weightsvector), 0, asInteger(maxNumDyadTypes)*sizeof(int));

  const char *outn[] = {"y", "x", "weightsvector", ""};
  SEXP outl = PROTECT(mkNamed(VECSXP, outn));
  SET_VECTOR_ELT(outl, 0, responsevec);
  SET_VECTOR_ELT(outl, 1, covmat);
  SET_VECTOR_ELT(outl, 2, weightsvector);

  MpleInit_hash_wl_RLE(s, INTEGER(responsevec), REAL(covmat), INTEGER(weightsvector), &wlm, asInteger(maxNumDyads), asInteger(maxNumDyadTypes));

  ErgmStateDestroy(s);
  PutRNGstate(); /* Must be called after GetRNGstate before returning to R */
  UNPROTECT(4);
  return outl;
}

/*************
Hashes the covariates, offset, and response onto an unsigned integer in the interval [0,numRows).
Uses Jenkins One-at-a-Time hash.

numRows should, ideally, be a power of 2, but doesn't have to be.
**************/
static inline unsigned int hashCovMatRow(double *newRow, unsigned int rowLength, unsigned int numRows,
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

void MpleInit_hash_wl_RLE(ErgmState *s, int *responsevec, double *covmat, int *weightsvector,
			  RLEBDM1D *wl, 
			  Edge maxNumDyads, Edge maxNumDyadTypes){
  Network *nwp = s->nwp;
  Model *m = s->m;

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
    
    for(Dyad i = 0; i < MIN(maxNumDyads,dc); i++, d=NextRLEBDM1D(d, step, wl, &r)){
      R_CheckUserInterruptEvery(1024u, i);
      Dyad2TH(&t, &h, d, N_NODES);
      
      int response = IS_OUTEDGE(t,h);
      ChangeStats1(t, h, nwp, m, response);
      if(response){
	for(unsigned int l=0; l<m->n_stats; l++){
	  m->workspace[l] = -m->workspace[l];
	}
      }
      /* In next command, if there is an offset vector then its total
	 number of entries should match the number of times through the 
	 inner loop (i.e., the number of dyads in the network) */          
      if(!insCovMatRow(m->workspace, covmat, m->n_stats,
		       maxNumDyadTypes, response, 
		       responsevec, weightsvector)) {
	warning("Too many unique dyads. MPLE is approximate, and MPLE standard errors are suspect.");
	break;
      }
    }
  } // End scope for loop variables.
}
