/*  File src/changestats_indices.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2019 Statnet Commons
 */
#include "changestats_indices.h"

/*****************                       
 changestat: d_indices
*****************/
C_CHANGESTAT_FN(c_indices) { 
  ZERO_ALL_CHANGESTATS(i);
    int edgeflag = IS_OUTEDGE(tail,head);
    CHANGE_STAT[0] += edgeflag ? -(int)tail : (int)tail;
    CHANGE_STAT[1] += edgeflag ? -(int)head : (int)head;
}
