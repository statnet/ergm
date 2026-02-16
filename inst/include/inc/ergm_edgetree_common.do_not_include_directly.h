/*  File inst/include/inc/ergm_edgetree_common.do_not_include_directly.h in
 *  package ergm, part of the Statnet suite of packages for network analysis,
 *  https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "../ergm_edgetree_types.h"
#include "../ergm_variadic_macros.h"

#define MIN(a,b) ((a)<(b) ? (a) : (b))
#define MAX(a,b) ((a)<(b) ? (b) : (a))

#define _DYADCOUNT1(nwp) _DYADCOUNT3(nwp->nnodes, nwp->bipartite, nwp->directed_flag)
#define _DYADCOUNT3(nnodes, bipartite, directed) ((bipartite)? (Dyad)((nnodes)-(bipartite))*(Dyad)(bipartite) : ((directed)? (Dyad)(nnodes)*(Dyad)((nnodes)-1) : (((Dyad)(nnodes)*(Dyad)((nnodes)-1))/2)))
#define DYADCOUNT(...) _GET_OVERRIDE3(__VA_ARGS__, _DYADCOUNT3, , _DYADCOUNT1,)(__VA_ARGS__)
#define EDGECOUNT(nwp) nwp->nedges

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
