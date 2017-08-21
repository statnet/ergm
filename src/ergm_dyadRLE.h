#ifndef _ERGM_DYADRLE_H_
#define _ERGM_DYADRLE_H_

#include "edgetree.h"

/* # Standard format for RLE-encoded dyad matrix for MHproposals #

   Let double *x point to the start of the dyad matrix information
   segment and let r be the number of runs of nonzero dyads. Here, the
   slice notation x[a : b] denotes elements of the array x[a] through
   x[b] inclusive.

   The matrix is required to be square. For bipartite and undirected
   networks, this simply means that only block- or triangle-subset of
   dyads are nonzero.

   x[0] = number of runs (r)

   x[1] = number of nonzero dyads (i.e., sum of all run lengths)

   x[2 : 2+r] = an r-vector whose 0th element is 1 and whose i'th
     element is the index of the first cell (listed in column-major
     order, indexed from 1) of a run.

   x[2+r+1 : 2+r+1+r+1] = an r+1-vector of cumulative run lengths, which also serve as scaled cumulative probabilities for
     selecting each run in the R code. (The 0th element
     is always 0, the last element is always
     sum of all run lengths, etc.)

*/

typedef unsigned int RLERun;

// A boolean (binary) square matrix encoded with run-length encoding, with run starts and lengths being represented as double-precision floats. 
typedef struct {
  Vertex n; // number of nodes
  RLERun nruns; // number of runs
  Dyad ndyads; // number of dyads
  double *starts; // start of a run of free dyads
  double *cumlens; // cumulative lengths of runs of free dyads
} BoolRLESqMatrixD;


/**
Unpack input from R into dyad matrix sampling information

Unpack a double *input vector containing block diagonal sampling
information into a `BoolRLESqMatrixD` structure and advance the
pointer to the end of the segment.

@param inputs a pointer to a pointer to a vector of inputs; will be
  updated by the procedure
@param network size
*/
static inline BoolRLESqMatrixD unpack_BoolRLESqMatrixD(double **inputs, Vertex n){
  double *x = *inputs;
  BoolRLESqMatrixD out = {
    .n = n,
    .ndyads = *(x++),
    .nruns = *(x++)
  };
  out.starts = x; x += out.nruns;
  out.cumlens = x; x += out.nruns+1;
  *inputs = x;
  return out;
}

/**
Generate a random dyad that belongs to a run

@param tail pointer to which to assign the tail value
@param head pointer to which to assign the head value
@param r RLE information
*/
static inline void GetRandDyadRLE(Vertex *tail, Vertex *head, const BoolRLESqMatrixD *m){
  // Select a dyad index at random
  Dyad d = unif_rand() * m->ndyads + 1;

  // Find the correct run via binary search.
  // FIXME: Radix search could be faster than O(log(nruns))?
  RLERun l=1, h=m->nruns;
  while(l!=h){
    RLERun r = (l+h)/2;
    if(d>m->cumlens[r]) l=r+1; else h=r;
  }

  // Dyad ID
  d = m->starts[l-1] + (d-m->cumlens[l-1]);

  // Quotion is the column (head) and remainder is the row (tail), counted from 0.
  ldiv_t q = ldiv(d-1, m->n);
  
  *tail = q.rem + 1;
  *head = q.quot + 1;
}


/**
Test if a dyad is present in a run

@param tail tail index (from 1)
@param head head index (from 1)
@param r RLE information
*/
static inline unsigned int GetDyadRLE(Vertex tail, Vertex head, const BoolRLESqMatrixD *m){
  Dyad d = (tail-1)+m->n*(Dyad)(head-1) + 1;

  // Find the correct run via binary search.
  // FIXME: Radix search could be faster than O(log(nruns))?
  RLERun l=1, h=m->nruns;
  while(l!=h){
    RLERun r = (l+h)/2;
    if(d<m->starts[r-1]) h=r-1; else l=r;
  }

  // If d's 0-indexed position in rth run is strictly less than the
  // length of that run, the dyad is nonzero.
  return (d-m->starts[l-1]) < m->cumlens[l]-m->cumlens[l-1];
}



#endif // _ERGM_DYADRLE_H_
