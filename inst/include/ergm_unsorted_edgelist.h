/*  File inst/include/ergm_unsorted_edgelist.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2022 Statnet Commons
 */
#ifndef _ERGM_UNSORTED_EDGELIST_H_
#define _ERGM_UNSORTED_EDGELIST_H_

#include "ergm_edgetree_types.h"
#include <R.h>

/*
   This is a data structure to store a collection of edges in unsorted order,
   using buffered arrays for tails and heads that offer amortized O(1) insertion
   (typically of a known-absent edge, to avoid duplicates), O(1) uniform random
   sampling, and O(1) deletion of the last sampled edge. It is often used as an
   auxiliary data structure to support stratified sampling, with one UnsrtEL for
   each stratum.

   For now, these are static-inlined for simplicity, but
   we may want to move some of the less frequently used functions
   to a C file and save RAM.
*/

typedef struct {
  Vertex *tails;
  Vertex *heads;
  unsigned int lindex;
  unsigned int nedges;
  unsigned int maxedges;
} UnsrtEL;

static inline UnsrtEL *UnsrtELInitialize(unsigned int nedges, Vertex *tails, Vertex *heads, Rboolean copy) {
  UnsrtEL *el = Calloc(1, UnsrtEL);

  if(nedges == 0) {
    el->tails = NULL;
    el->heads = NULL;
  } else {
    if(copy) {
      el->tails = Calloc(nedges, Vertex);
      memcpy(el->tails, tails, nedges*sizeof(Vertex));
      el->heads = Calloc(nedges, Vertex);
      memcpy(el->heads, heads, nedges*sizeof(Vertex));
    } else {
      el->tails = tails;
      el->heads = heads;
    }
    el->tails--; // Index from 1.
    el->heads--; // Index from 1.
  }
  el->lindex = 0;
  el->nedges = nedges;
  el->maxedges = nedges;

  return el;
}

static inline void UnsrtELDestroy(UnsrtEL *el) {
  if(el->tails) {
    el->tails++;
    el->heads++;
    Free(el->tails);
    Free(el->heads);
  }
  Free(el);
}

static inline void UnsrtELClear(UnsrtEL *el) {
  el->lindex = 0;
  el->nedges = 0;
}

static inline unsigned int UnsrtELGetRand(Vertex *tail, Vertex *head, UnsrtEL *el) {
  if(el->nedges == 0) {
    return 0;
  }

  el->lindex = el->nedges*unif_rand() + 1;
  *tail = el->tails[el->lindex];
  *head = el->heads[el->lindex];
  return el->lindex;
}

static inline unsigned int UnsrtELSearch(Vertex tail, Vertex head, UnsrtEL *el) {
  if(el->lindex > 0 && tail == el->tails[el->lindex] && head == el->heads[el->lindex]) {
   // found edge at last index
   return el->lindex;
  } else {
    // linear search

#ifdef DEBUG_UnsrtEL
    /* To stop on an inefficient search, have the debugger break on the next line. */
    Rprintf("UnsrtELSearch() called for an edge other than the last one inserted or selected, resulting in a linear search. This is O(E) slow and should be avoided whenever possible.\n");
#endif

    unsigned int i = el->nedges;
    while(i != 0 && (tail != el->tails[i] || head != el->heads[i])) {
      i--;
    }

    if(i) {
      el->lindex = i;
    }
    return i;
  }
}

static inline void UnsrtELDelete(Vertex tail, Vertex head, UnsrtEL *el) {
  if(UnsrtELSearch(tail, head, el) == 0) {
    return; // edge already absent
  }

  el->tails[el->lindex] = el->tails[el->nedges];
  el->heads[el->lindex] = el->heads[el->nedges];
  el->nedges--;
  el->lindex = 0;
}

// used by HashEL; assumes `at` is a valid edge index in `el`
static inline void UnsrtELDeleteAt(unsigned int at, UnsrtEL *el) {
  if(at < el->nedges) {
    el->tails[at] = el->tails[el->nedges];
    el->heads[at] = el->heads[el->nedges];
  }

  el->nedges--;
  el->lindex = 0;
}

static inline void UnsrtELInsert(Vertex tail, Vertex head, UnsrtEL *el) {
  if(el->nedges == el->maxedges) {
    el->maxedges = MAX(el->maxedges, 1) * 2;
    if(el->tails) {
      el->tails++;
      el->heads++;
    }
    el->tails = Realloc(el->tails, el->maxedges, Vertex) - 1;
    el->heads = Realloc(el->heads, el->maxedges, Vertex) - 1;
  }

#ifdef DEBUG_UnsrtEL
  unsigned int i = el->nedges;
  while(i != 0 && (tail != el->tails[i] || head != el->heads[i])) {
    i--;
  }
  if(i) { /* To stop on an invariant violation, have the debugger break on the next line. */
    Rprintf("UnsrtELInsert() called for an edge already present in the list. This should never happen.\n");
  }
#endif

  el->nedges++;
  el->lindex = el->nedges;
  el->tails[el->lindex] = tail;
  el->heads[el->lindex] = head;
}

// helper function for cleaner client code
static inline void UnsrtELToggleKnown(Vertex tail, Vertex head, UnsrtEL *el, int edgestate) {
  if(edgestate) {
    UnsrtELDelete(tail, head, el);
  } else {
    UnsrtELInsert(tail, head, el);
  }
}

#endif // _ERGM_UNSORTED_EDGELIST_H_
