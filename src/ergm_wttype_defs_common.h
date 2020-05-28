#define PROP_PRINT Rprintf("  (%d, %d) -> %f  ", MHp->toggletail[i], MHp->togglehead[i], MHp->toggleweight[i])
#define PROP_CHANGESTATS WtChangeStats(MHp->ntoggles, MHp->toggletail, MHp->togglehead, MHp->toggleweight, nwp, m)
#define PROP_COMMIT WtSetEdge(MHp->toggletail[i], MHp->togglehead[i], MHp->toggleweight[i], nwp)
#define DISPATCH_ErgmState WtErgmState
#define DISPATCH_ErgmStateInit WtErgmStateInit
#define DISPATCH_Model WtModel
#define DISPATCH_MHProposal WtMHProposal
#define DISPATCH_ErgmStateRSave WtErgmStateRSave
#define DISPATCH_ErgmStateDestroy WtErgmStateDestroy
#define DISPATCH_MCMCSample WtMCMCSample
#define DISPATCH_Network WtNetwork
