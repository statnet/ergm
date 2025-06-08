/*  File inst/include/ergm_wtedgetree.h in package ergm, part of the Statnet
 *  suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#ifndef _ERGM_WTEDGETREE_H_
#define _ERGM_WTEDGETREE_H_

#include "ergm_edgetype_set_double.h"

#include "inc/ergm_edgetree.h.template.do_not_include_directly.h"

#include "ergm_edgetype_unset.h"

/*
  Workaround to enable WtNetworkInitialize() to be called with either
  7 or 10 arguments.

  TODO: Remove the workaround and rename WtNetworkInitialize_noLT() to
    WtNetworkInitialize() after no CRAN packages use it.
*/

#define WtNetworkInitialize(...) _GET_OVERRIDE10(__VA_ARGS__, WtNetworkInitialize_LT, , , WtNetworkInitialize_noLT, , , , , , , , )(__VA_ARGS__)

#define WtNetworkInitialize_LT(tails, heads, weights, nedges, nnodes, directed_flag, bipartite, lasttoggle_flag, time, lasttoggle) \
  WtNetworkInitialize_noLT((tails), (heads), (weights), (nedges), (nnodes), (directed_flag), (bipartite))

#endif
