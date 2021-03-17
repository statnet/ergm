#define PROP_PRINT Rprintf("  (%d, %d)  ", MHp->toggletail[i], MHp->togglehead[i])
#define PROP_CHANGESTATS ChangeStats(MHp->ntoggles, MHp->toggletail, MHp->togglehead, nwp, m)
#define PROP_COMMIT ToggleEdge(MHp->toggletail[i], MHp->togglehead[i], nwp)
#define DISPATCH_ErgmState ErgmState
#define DISPATCH_ErgmStateInit ErgmStateInit
#define DISPATCH_Model Model
#define DISPATCH_MHProposal MHProposal
#define DISPATCH_ErgmStateRSave ErgmStateRSave
#define DISPATCH_ErgmStateDestroy ErgmStateDestroy
#define DISPATCH_Network Network
