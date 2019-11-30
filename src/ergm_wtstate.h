#ifndef _ERGM_WTSTATE_H_
#define _ERGM_WTSTATE_H_

#include "ergm_wtedgetree.h"
#include "ergm_wtchangestat.h"
#include "ergm_wtMHproposal.h"
#include "ergm_wtmodel.h"

typedef struct{
  WtNetwork *nwp;
  WtModel *m;
  WtMHProposal *MHp;
} ErgmWtState;

ErgmWtState *ErgmWtStateInit(// Network settings
			     Vertex n_nodes, Rboolean directed_flag, Vertex bip,
                             // Model settings
			     int nterms, const char *funnames, const char *sonames,
                             // Proposal settings
			     const char *MHProposaltype, const char *MHProposalpackage,
			     double *inputs,
                             // Network state
                             Edge n_edges,
			     Vertex *tails, Vertex *heads, double *weights);
void ErgmWtStateDestroy(ErgmWtState *s);

#endif
