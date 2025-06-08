/*  File src/changestats_dyad_ind.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#include "ergm_changestat.h"
#include "ergm_edgetype_set_binary.h"
#define ECHANGE(a) (edgestate ? -(a) : +(a))
#define ECHANGE1 (edgestate ? -1 : +1)
#define SVARIANT(a) a
#include "changestats_dyad_ind.c.template.do_not_include_directly.h"

S_CHANGESTAT_FN(s_edges) {
  CHANGE_STAT[0] = N_EDGES;
}
