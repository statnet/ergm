/*  File src/ergm_type_defs_common.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */

#define EDGETYPE(FUN) FUN

#define PROP_PRINT Rprintf("  (%d, %d)  ", MHp->toggletail[i], MHp->togglehead[i])
#define PROP_CHANGESTATS ChangeStats(MHp->ntoggles, MHp->toggletail, MHp->togglehead, nwp, m)
#define PROP_CHANGESTATS_DO ChangeStatsDo(MHp->ntoggles, MHp->toggletail, MHp->togglehead, nwp, m)
#define PROP_CHANGESTATS_UNDO ChangeStatsUndo(MHp->ntoggles, MHp->toggletail, MHp->togglehead, nwp, m)
#define PROP_COMMIT ToggleEdge(MHp->toggletail[i], MHp->togglehead[i], nwp)
#define PROP_FINISH ToggleEdge(MHp->toggletail[MHp->ntoggles-1], MHp->togglehead[MHp->ntoggles-1], nwp)
