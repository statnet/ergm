#ifndef _ERGM_UNSORTED_EDGELIST_H_
#define _ERGM_UNSORTED_EDGELIST_H_

#include "ergm_edgetree_types.h"
#include <R.h>

/* This is a data structure based on the StratTNT algorithm developed
   by Chad Klumb. For now, these are static-inlined for simplicity, but
   we may want to move some of the less frequently used functions
   to a C file and save RAM.*/

typedef struct{
  Vertex *tails, *heads;
  Edge lindex, nedges, maxedges;
} UnsrtEL;


static inline UnsrtEL *UnsrtELInitialize(Edge nedges, Vertex *tails, Vertex *heads, Rboolean copy){
  UnsrtEL *el = Calloc(1, UnsrtEL);

  Vertex *mytails, *myheads;
  if(nedges==0) mytails = myheads = NULL; // Realloc will handle this correctly on insertion.
  else if(copy){
    mytails = Calloc(nedges, Vertex);
    memcpy(mytails, tails, nedges*sizeof(Vertex));
    myheads = Calloc(nedges, Vertex);
    memcpy(myheads, heads, nedges*sizeof(Vertex));
  }else{
    mytails = tails;
    myheads = heads;
  }
  el->tails = mytails-1; // Index from 1.
  el->heads = myheads-1; // Index from 1.
  el->lindex = 0;
  el->nedges = el->maxedges = nedges;

  return el;
}

static inline void UnsrtELDestroy(UnsrtEL *el){
  el->tails++; Free(el->tails);
  el->heads++; Free(el->heads);
  Free(el);
}

static inline void UnsrtELClear(UnsrtEL *el) {
  el->lindex = 0;
  el->nedges = 0;
}

static inline Edge UnsrtELGetRand(Vertex *tail, Vertex *head, UnsrtEL *el){
  if(el->nedges==0) return 0;
  el->lindex = el->nedges*unif_rand() + 1;
  *tail = el->tails[el->lindex];
  *head = el->heads[el->lindex];
  return el->lindex;
}

static inline Edge UnsrtELSearch(Vertex tail, Vertex head, UnsrtEL *el){
  if(el->lindex==0 || tail!=el->tails[el->lindex] || head!=el->heads[el->lindex]){ // Linear search
#ifdef DEBUG_UnsrtEL
    /* To stop on an inefficient search, have the debugger break on the next line. */
    Rprintf("UnsrtELSearch() called for an edge other than the last one inserted or selected, resulting in a linear search. This is O(E) slow and should be avoided whenever possible.\n");
#endif
    Edge i = el->nedges;
    while(i!=0 && (tail!=el->tails[i] || head!=el->heads[i])) i--;

    if(i) return el->lindex=i;
    else return 0;
  }

  return el->lindex;
}

static inline void UnsrtELDelete(Vertex tail, Vertex head, UnsrtEL *el){
  if(UnsrtELSearch(tail, head, el)==0) return; // Edge already absent.
  el->tails[el->lindex] = el->tails[el->nedges];
  el->heads[el->lindex] = el->heads[el->nedges];
  el->nedges--;
  el->lindex=0;
}

// used by HashEL; assumes `at` is a valid edge index in `el`
static inline void UnsrtELDeleteAt(unsigned int at, UnsrtEL *el) {
  if(at < el->nedges) {
    el->tails[at] = el->tails[el->nedges];
    el->heads[at] = el->heads[el->nedges];
  }
  
  el->nedges--;
  el->lindex=0;
}

static inline void UnsrtELInsert(Vertex tail, Vertex head, UnsrtEL *el){
  if(el->nedges==el->maxedges){
    el->maxedges = MAX(el->maxedges,1) * 2;
    el->tails = Realloc(el->tails+1, el->maxedges, Vertex) - 1;
    el->heads = Realloc(el->heads+1, el->maxedges, Vertex) - 1;
  }

#ifdef DEBUG_UnsrtEL
  Edge i = el->nedges;
  while(i!=0 && (tail!=el->tails[i] || head!=el->heads[i])) i--;
  if(i) /* To stop on an invariant violation, have the debugger break on the next line. */
    Rprintf("UnsrtELInsert() called for an edge already present in the list. This should never happen.\n");
#endif

  el->lindex = ++el->nedges;
  el->tails[el->lindex] = tail;
  el->heads[el->lindex] = head;
}

// helper function for cleaner client code
static inline void UnsrtELToggleKnown(Vertex tail, Vertex head, UnsrtEL *el, int edgeflag) {
  if(edgeflag) {
    UnsrtELDelete(tail, head, el);
  } else {
    UnsrtELInsert(tail, head, el);
  }
}

#endif // _ERGM_UNSORTED_EDGELIST_H_
