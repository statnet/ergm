/*  File src/edgetree.c in package ergm, part of the Statnet suite of packages
 *  for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#include "ergm_edgetree.h"
#include "ergm_edgetype_set_binary.h"
#include "edgetree.c.template.do_not_include_directly.h"

/*****************
 Edge ToggleKnownEdge

 Toggle an edge whose status is known:  Set it to the opposite of its current
 value.  Return 1 if edge added, 0 if deleted.
*****************/

/* *** don't forget tail->head, so this function now accepts tail before head */

void ToggleKnownEdge (Vertex tail, Vertex head, Network *nwp, Rboolean edgestate)
{
#ifndef NDEBUG
  if((EdgetreeSearch(tail, head, nwp->outedges)!=0) != edgestate)
    error("ToggleKnownEdge() called with an incorrect edgestate. Note that this produces an error only if compiling with NDEBUG macro unset and silently produces undefined behavior otherwise.");
#endif // NDEBUG
  ENSURE_TH_ORDER;
  if (edgestate){
    DeleteEdgeFromTrees(tail,head,nwp);
  }else{
    AddEdgeToTrees(tail,head,nwp);
  }
}
