/*  File src/ergm_wttype_defs_common.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#define PROP_PRINT Rprintf("  (%d, %d) -> %f  ", MHp->toggletail[i], MHp->togglehead[i], MHp->toggleweight[i])
#define PROP_CHANGESTATS WtChangeStats(MHp->ntoggles, MHp->toggletail, MHp->togglehead, MHp->toggleweight, nwp, m)
#define PROP_CHANGESTATS_DO WtChangeStatsDo(MHp->ntoggles, MHp->toggletail, MHp->togglehead, MHp->toggleweight, nwp, m)
#define PROP_CHANGESTATS_UNDO WtChangeStatsUndo(MHp->ntoggles, MHp->toggletail, MHp->togglehead, MHp->toggleweight, nwp, m)
#define PROP_COMMIT WtSetEdge(MHp->toggletail[i], MHp->togglehead[i], MHp->toggleweight[i], nwp)
#define PROP_FINISH WtSetEdge(MHp->toggletail[MHp->ntoggles-1], MHp->togglehead[MHp->ntoggles-1], MHp->toggleweight[MHp->ntoggles-1], nwp)
#define EDGETYPE_ErgmState WtErgmState
#define EDGETYPE_ErgmStateInit WtErgmStateInit
#define EDGETYPE_Model WtModel
#define EDGETYPE_MHProposal WtMHProposal
#define EDGETYPE_ErgmStateRSave WtErgmStateRSave
#define EDGETYPE_ErgmStateDestroy WtErgmStateDestroy
#define EDGETYPE_Network WtNetwork
