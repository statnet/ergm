/*  File inst/include/ergm_rlebdm.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#ifndef _ERGM_RLEBDM_H_
#define _ERGM_RLEBDM_H_

#include <limits.h>
#include "ergm_edgetree.h"

/* Serialization format for RLE-encoded Binary Dyad Matrix with only
   TRUE (1) values stored, with indices and lengths stored as Double
   (RLEBDM1D), serialized as a double type vector.

   Let double *x point to the start of the dyad matrix information
   segment and let r be the number of runs of TRUE dyads. Here, the
   slice notation x[a : b] denotes elements of the array x[a] through
   x[b] inclusive.

   The matrix is required to be square. For bipartite and undirected
   networks, this simply means that only block- or triangle-subset of
   dyads are TRUE.

   x[0] = number of nodes in the network

   x[1] = number of TRUE dyads (i.e., sum of all run lengths)

   x[2 : 2+r] = an r-vector whose i'th
     element is the index of the first cell (listed in column-major
     order, indexed from 1) of the i+1st run.

   x[2+r+1 : 2+r+1+r+1] = an r+1-vector of cumulative run lengths, which also serve as scaled cumulative probabilities for
     selecting each run in the R code. (The 0th element
     is always 0, the last element is always
     sum of all run lengths, etc.)

*/

#define TH2Dyad(n, tail, head) ((tail-1)+n*(Dyad)(head-1) + 1)
#define Dyad2T(n, d) ((Vertex)((d-1)%(Vertex)n) + 1)
#define Dyad2H(n, d) ((Vertex)(((Dyad)(d-1))/(Vertex)n) + 1)

static inline void Dyad2TH(Vertex *tail, Vertex *head, Dyad d, Vertex n){
  // Quotient is the column (head) and remainder is the row (tail), counted from 0.
  // On systems for which long integers are 32 bits, use lldiv() instead of ldiv().
#if ULONG_MAX >= 18446744073709551615ULL
  ldiv_t q = ldiv(d-1, n);
#else
  lldiv_t q = lldiv(d-1, n);
#endif
  *tail = q.rem + 1;
  *head = q.quot + 1;
}

typedef unsigned int RLERun;

// A boolean (binary) square matrix encoded with run-length encoding, with run starts and lengths being represented as double-precision floats. 
typedef struct {
  Vertex n; // number of nodes
  RLERun nruns; // number of runs
  Dyad ndyads; // number of dyads
  double *starts; // start of a run of free dyads
  double *cumlens; // cumulative lengths of runs of free dyads
  unsigned int maxlen; // length of the longest run of free dyads (useful in rejection sampling)
} RLEBDM1D;


/**
Unpack input from R into dyad matrix sampling information

Unpack a double *input vector containing block diagonal sampling
information into a `RLEBDM1D` structure and advance the
pointer to the end of the segment.

@param inputs a pointer to a pointer to a vector of inputs; will be
  updated by the procedure
@param network size
*/
static inline RLEBDM1D unpack_RLEBDM1D(double **inputs){
  double *x = *inputs;
  RLEBDM1D out;
  out.n = *(x++);
  out.ndyads = *(x++);
  out.nruns = *(x++);
  out.starts = x; x += out.nruns;
  out.cumlens = x; x += out.nruns+1;
  *inputs = x;

  out.maxlen=0;
  for(unsigned int r=1; r<=out.nruns; r++){
    unsigned int l = out.cumlens[r]-out.cumlens[r-1];
    if(l > out.maxlen) out.maxlen = l;
  }
  
  return out;
}

/**
Generate a random dyad that belongs to a run using inverse transform sampling

@param tail pointer to which to assign the tail value
@param head pointer to which to assign the head value
@param r RLE information

@note This procedure is O(log(n)) in m->nruns.
*/
static inline void GetRandRLEBDM1D_ITS(Vertex *tail, Vertex *head, const RLEBDM1D *m){
  // Select a dyad index at random
  Dyad i = unif_rand() * m->ndyads + 1;

  // Find the correct run via binary search.
  // FIXME: Radix search could be faster than O(log(nruns))?
  RLERun l=1, h=m->nruns;
  while(l!=h){
    RLERun r = (l+h)/2;
    if(i>m->cumlens[r]) l=r+1; else h=r;
  }

  unsigned int rl = m->cumlens[l] - m->cumlens[l-1];
  
  // Dyad ID
  Dyad d = (Dyad)m->starts[l-1] + (rl==1 ? 0 : rl*unif_rand());

  Dyad2TH(tail, head, d, m->n);
}


/* Set the default generator. */
#define GetRandRLEBDM1D GetRandRLEBDM1D_RS

/**
Generate a random dyad that belongs to a run using rejection sampling

@param tail pointer to which to assign the tail value
@param head pointer to which to assign the head value
@param r RLE information

@note This procedure is O(1) in m->nruns, but with a larger constant
due to multiple unif_rand() calls. However, because runs tend to have
similar lengths, the acceptance probability is likely to be high.

*/
static inline void GetRandRLEBDM1D_RS(Vertex *tail, Vertex *head, const RLEBDM1D *m){
  RLERun r;
  double l;
  double u;
  double pr;
  do{
    // Select a run at random
    u = unif_rand(); // ~Uniform(0,1)
    double x = u * m->nruns + 1; // ~Uniform(1, m->nruns+1)
    r = floor(x); // ~DUniform(1..m->nruns)
    u = x - r; // ~Uniform(0,1)
    l = m->cumlens[r]-m->cumlens[r-1];
    pr = l/m->maxlen;
  }while(pr < u);
  // So u < pr.

  // u /= pr; // ~Uniform(0,1); but there might not be enough bits left.

  // The run r now has probability of selection proportional to its
  // length, and l is its length.

  // Dyad ID
  Dyad d = (Dyad)m->starts[r-1] + (l==1 ? 0 : unif_rand()*l);

  Dyad2TH(tail, head, d, m->n);
}

/**
Test if a dyad is present in a run

@param tail tail index (from 1)
@param head head index (from 1)
@param m RLE information
*/
static inline unsigned int GetRLEBDM1D(Vertex tail, Vertex head, const RLEBDM1D *m){
  Dyad d = TH2Dyad(m->n, tail, head);

  if(d<m->starts[0]) return FALSE; // d precedes the first run.
  
  // Find the correct run via binary search.
  // FIXME: Radix search could be faster than O(log(nruns))?
  RLERun l=1, h=m->nruns;
  while(l!=h){
    RLERun r = (l+h+1)/2;
    if(d<m->starts[r-1]) h=r-1; else l=r;
  }

  // If d's 0-indexed position in rth run is strictly less than the
  // length of that run, the dyad is TRUE.
  return (d-m->starts[l-1]) < m->cumlens[l]-m->cumlens[l-1];
}

/**
Advance to the strideth next TRUE dyad index 

@param d current dyad's number
@param stride how much to advance
@param m dyad matrix
@param hint a pointer to an integer containing the RLE run, which is updated; set to 0 to initialise and to NULL to disable
*/
static inline Dyad NextRLEBDM1D(Dyad d, Dyad stride, const RLEBDM1D *m, RLERun *hint){
  RLERun l;
  if(hint && *hint!=0){
    l = *hint;
  }else{
    l = 1;
    RLERun h=m->nruns;
    while(l!=h){
      RLERun r = (l+h+1)/2;
      if(d<m->starts[r-1]) h=r-1; else l=r;
    }
  }

  // 0-based index of the current dyad in the list of TRUE dyads.
  Dyad di = d - m->starts[l-1] + m->cumlens[l-1];
  // 0-based index of the next dyad in the list of TRUE dyads.
  Dyad nxtdi = (di+stride) % m->ndyads;

  RLERun nxtl;  
  if(nxtdi < di){ // We've wrapped around.
    nxtl = 1;
  }else{ // We haven't wrapped around.
    nxtl = l;
  }

  while(nxtdi >= m->cumlens[nxtl]) nxtl++;
  // nxtl is now the 1-based run index of the next dyad.

  // Update the hint if active.
  if(hint) *hint = nxtl;
  
  return nxtdi - m->cumlens[nxtl-1] + m->starts[nxtl-1];
}

#define FirstRLEBDM1D(m) ((m)->starts[0])
#define LastRLEBDM1D(m) ((m)->starts[(m)->nruns-1] - (m)->cumlens[(m)->nruns-2] + (m)->cumlens[(m)->nruns-1] - 1)

void PrintRLEBDM1D(const RLEBDM1D *m);

#endif // _ERGM_RLEBDM_H_
