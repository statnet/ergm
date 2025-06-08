/*  File inst/include/ergm_edgetree.h in package ergm, part of the Statnet
 *  suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#ifndef _ERGM_EDGETREE_H_
#define _ERGM_EDGETREE_H_

#include "ergm_edgetype_set_binary.h"

#include "inc/ergm_edgetree.h.template.do_not_include_directly.h"

#include "ergm_edgetype_unset.h"

void ToggleKnownEdge (Vertex tail, Vertex head, Network *nwp, Rboolean edgestate);

/*
  Workaround to enable NetworkInitialize() to be called with either 6
  or 9 arguments.

  TODO: Remove the workaround and rename NetworkInitialize_noLT() to
    NetworkInitialize() after no CRAN packages use it.
*/
#define NetworkInitialize(...) _GET_OVERRIDE9(__VA_ARGS__, NetworkInitialize_LT, , , NetworkInitialize_noLT, , , , , , , , )(__VA_ARGS__)

#define NetworkInitialize_LT(tails, heads, nedges, nnodes, directed_flag, bipartite, lasttoggle_flag, time, lasttoggle) \
  NetworkInitialize_noLT((tails), (heads), (nedges), (nnodes), (directed_flag), (bipartite))

#endif
