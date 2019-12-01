#ifndef _ERGM_WTSTATE_H_
#define _ERGM_WTSTATE_H_

#include "ergm_wtedgetree.h"
#include "ergm_wtchangestat.h"
#include "ergm_wtMHproposal.h"
#include "ergm_wtmodel.h"

#define NO_WTMODEL 0, NULL, NULL
#define NO_WTMHPROPOSAL  NULL, NULL
#define NO_LASTTOGGLE FALSE, 0, NULL
#define NO_WTNWSTATE 0, NULL, NULL, NULL

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
			     Vertex *tails, Vertex *heads, double *weights,
                             Rboolean timings, int time, int *lasttoggle);
void ErgmWtStateDestroy(ErgmWtState *s);

#endif
