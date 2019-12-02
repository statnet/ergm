#ifndef _ERGM_STATE_H_
#define _ERGM_STATE_H_

#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_MHproposal.h"
#include "ergm_model.h"

#define NO_MODEL 0, NULL, NULL, FALSE
#define NO_MHPROPOSAL  NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0
#define NO_LASTTOGGLE FALSE, 0, NULL
#define NO_NWSTATE 0, NULL, NULL

typedef struct{
  Network *nwp;
  Model *m;
  MHProposal *MHp;
} ErgmState;

ErgmState *ErgmStateInit(// Network settings
                         Vertex n_nodes, Rboolean directed_flag, Vertex bip,
                         // Model settings
                         int nterms, const char *funnames, const char *sonames, Rboolean noinit_s,
                         // Proposal settings
                         const char *MHProposaltype, const char *MHProposalpackage,
                         int *attribs, int *maxout, int *maxin, int *minout,
                         int *minin, int condAllDegExact, int attriblength,
                         // Numeric inputs
                         double *inputs,
                         // Network state
                         Edge n_edges,
                         Vertex *tails, Vertex *heads,
                         Rboolean timings, int time, int *lasttoggle);
void ErgmStateDestroy(ErgmState *s);

#endif
