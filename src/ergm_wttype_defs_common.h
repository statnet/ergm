/*  File src/ergm_wttype_defs_common.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2023 Statnet Commons
 */
#define PROP_PRINT Rprintf("  (%d, %d) -> %f  ", MHp->toggletail[i], MHp->togglehead[i], MHp->toggleweight[i])
#define PROP_CHANGESTATS WtChangeStats(MHp->ntoggles, MHp->toggletail, MHp->togglehead, MHp->toggleweight, nwp, m)
#define PROP_COMMIT WtSetEdge(MHp->toggletail[i], MHp->togglehead[i], MHp->toggleweight[i], nwp)
#define DISPATCH_ErgmState WtErgmState
#define DISPATCH_ErgmStateInit WtErgmStateInit
#define DISPATCH_Model WtModel
#define DISPATCH_MHProposal WtMHProposal
#define DISPATCH_ErgmStateRSave WtErgmStateRSave
#define DISPATCH_ErgmStateDestroy WtErgmStateDestroy
#define DISPATCH_Network WtNetwork
