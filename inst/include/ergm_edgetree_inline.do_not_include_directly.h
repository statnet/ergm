/*  File inst/include/ergm_edgetree_inline.do_not_include_directly.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */

#define EDGETYPE(FUN) FUN

#include "ergm_edgetree_inline_template.do_not_include_directly.h"

#undef EDGETYPE

/*****************
 int GetEdge

Get edge value. Return 0 if edge does not exist.
*****************/
static inline unsigned int GetEdge (Vertex tail, Vertex head, Network *nwp) 
{
  ENSURE_TH_ORDER;

  return EdgetreeSearch(tail,head,nwp->outedges)!=0;
}

