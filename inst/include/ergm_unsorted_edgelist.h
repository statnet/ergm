#ifndef _ERGM_UNSORTED_EDGELIST_H_
#define _ERGM_UNSORTED_EDGELIST_H_

#include "ergm_edgetree_types.h"
#include <R.h>

/* This is a data structure based on the StratTNT algorithm developed
   by Chad Kumb. For now, these are static-inlined for simplicity, but
   we may want to move some of the less frequently used functions
   to a C file and save RAM.*/

typedef struct{
  Vertex *tails, *heads, ltail, lhead;
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

static inline Edge UnsrtELGetRand(Vertex *tail, Vertex *head, UnsrtEL *el){
  if(el->nedges==0) return 0;
  el->lindex = el->nedges*unif_rand() + 1;
  *tail = el->ltail = el->tails[el->lindex];
  *head = el->lhead = el->heads[el->lindex];
  return el->lindex;
}

static inline Edge UnsrtELSearch(Vertex tail, Vertex head, UnsrtEL *el){
  if(el->lindex==0 || el->lindex>el->nedges || tail!=el->ltail || head!=el->lhead){ // Linear search
    Edge i = el->nedges;
    while(i!=0 && (tail!=el->tails[i] || head!=el->heads[i])) i--;

    if(i){
      el->ltail=tail;
      el->lhead=head;
      el->lindex=i;
    }else return 0;
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

static inline void UnsrtELInsert(Vertex tail, Vertex head, UnsrtEL *el){
  if(el->nedges==el->maxedges){
    el->maxedges = MAX(el->maxedges,1) * 2;
    el->tails = Realloc(el->tails+1, el->maxedges, Vertex) - 1;
    el->heads = Realloc(el->heads+1, el->maxedges, Vertex) - 1;
  }

  el->lindex = ++el->nedges;
  el->tails[el->lindex] = tail;
  el->heads[el->lindex] = head;
}

#endif // _ERGM_UNSORTED_EDGELIST_H_
