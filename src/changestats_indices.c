/*  File src/changestats_indices.c in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2021 Statnet Commons
 */
#include "ergm_edgetree.h"
#include "ergm_changestat.h"

/*****************                       
 changestat: d_indices
*****************/
C_CHANGESTAT_FN(c_indices) { 
    CHANGE_STAT[0] += edgestate ? -(int)tail : (int)tail;
    CHANGE_STAT[1] += edgestate ? -(int)head : (int)head;
}
