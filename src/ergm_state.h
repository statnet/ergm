#ifndef _ERGM_STATE_H_
#define _ERGM_STATE_H_

#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_MHproposal.h"
#include "ergm_model.h"

typedef struct{
  Network *nwp;
  Model *m;
  MHProposal *MHp;
} ErgmState;

ErgmState *ErgmStateInit(// Network settings
                         Vertex n_nodes, Rboolean directed_flag, Vertex bip,
                         // Model settings
                         int nterms, const char *funnames, const char *sonames,
                         // Proposal settings
                         const char *MHProposaltype, const char *MHProposalpackage,
                         int *attribs, int *maxout, int *maxin, int *minout,
                         int *minin, int condAllDegExact, int attriblength,
                         // Numeric inputs
                         double *inputs,
                         // Network state
                         Edge n_edges,
                         Vertex *tails, Vertex *heads);
void ErgmStateDestroy(ErgmState *s);

#endif
