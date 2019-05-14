/*  File inst/include/ergm_edgetree_common.do_not_include_directly.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2019 Statnet Commons
 */
#ifndef _ERGM_EDGETREE_COMMON_H_
#define  _ERGM_EDGETREE_COMMON_H_
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#define MIN(a,b) ((a)<(b) ? (a) : (b))
#define MAX(a,b) ((a)<(b) ? (b) : (a))

#define _GET_OVERRIDE1(_1, NAME, ...) NAME
#define _GET_OVERRIDE2(_1, _2, NAME, ...) NAME
#define _GET_OVERRIDE3(_1, _2, _3, NAME, ...) NAME
#define _GET_OVERRIDE4(_1, _2, _3, _4, NAME, ...) NAME
#define _GET_OVERRIDE5(_1, _2, _3, _4, _5, NAME, ...) NAME
#define _GET_OVERRIDE6(_1, _2, _3, _4, _5, _6, NAME, ...) NAME

#define _DYADCOUNT1(nwp) _DYADCOUNT3(nwp->nnodes, nwp->bipartite, nwp->directed_flag)
#define _DYADCOUNT3(nnodes, bipartite, directed) ((bipartite)? (unsigned long)((nnodes)-(bipartite))*(unsigned long)(bipartite) : ((directed)? (unsigned long)(nnodes)*(unsigned long)((nnodes)-1) : (((unsigned long)(nnodes)*(unsigned long)((nnodes)-1))/2)))
#define DYADCOUNT(...) _GET_OVERRIDE3(__VA_ARGS__, _DYADCOUNT3, , _DYADCOUNT1,)(__VA_ARGS__)
#define EDGECOUNT(nwp) nwp->nedges

/*typedef unsigned int Vertex; */
typedef unsigned int Vertex;
typedef unsigned int Edge;
typedef unsigned long int Dyad;

/* Ensure that tail < head for undriected networks. */
#define ENSURE_TH_ORDER							\
  if(!(nwp->directed_flag) && tail>head){				\
    Vertex temp;							\
    temp = tail;							\
    tail = head;							\
    head = temp;							\
  }

// This one is implemented as a macro, since it's very simple and works exactly the same for weighted and unweighted.
#define GetRandDyad(tail, head, nwp)					\
  if((nwp)->bipartite){							\
    *(tail) = 1 + unif_rand() * (nwp)->bipartite;			\
    *(head) = 1 + (nwp)->bipartite + unif_rand() * ((nwp)->nnodes - (nwp)->bipartite); \
  }else{								\
    *(tail) = 1 + unif_rand() * (nwp)->nnodes;				\
    *(head) = 1 + unif_rand() * ((nwp)->nnodes-1);			\
    if(*(head)>=*(tail)) (*(head))++;					\
    									\
    if (!(nwp)->directed_flag && *(tail) > *(head)) {			\
      Vertex tmp = *(tail);						\
      *(tail) = *(head);						\
      *(head) = tmp;							\
    }									\
  }

/* Dur_Inf is a structure containing information about durations of
edges in a network structure.
*/ 
typedef struct Dur_Infstruct {
  int time;
  int *lasttoggle;
} Dur_Inf;

#endif // _ERGM_EDGETREE_COMMON_H_
