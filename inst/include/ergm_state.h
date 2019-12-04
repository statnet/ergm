#ifndef _ERGM_STATE_H_
#define _ERGM_STATE_H_

#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_MHproposal.h"
#include "ergm_model.h"

#define ARGS_NWSETTINGS SEXP dn, SEXP dflag, SEXP bipartite
#define ARGS_MODEL SEXP mR
#define ARGS_MHPROPOSAL SEXP pR
#define ARGS_NWSTATE SEXP nedges, SEXP tails, SEXP heads
#define ARGS_LASTTOGGLE SEXP time, SEXP lasttoggle

#define YES_NWSETTINGS asInteger(dn), asInteger(dflag), asInteger(bipartite)
#define YES_MODEL mR, FALSE
#define YES_MODEL_NOINIT_S mR, TRUE
#define YES_MHPROPOSAL pR
#define YES_NWSTATE asInteger(nedges), (Vertex*) INTEGER(tails), (Vertex*) INTEGER(heads)
#define YES_LASTTOGGLE length(time)>0, asInteger(time), INTEGER(lasttoggle)

#define NO_MODEL NULL, FALSE
#define NO_MHPROPOSAL  NULL
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
                         SEXP mR, Rboolean noinit_s,
                         // Proposal settings
                         SEXP pR,
                         // Network state
                         Edge n_edges,
                         Vertex *tails, Vertex *heads,
                         Rboolean timings, int time, int *lasttoggle);
void ErgmStateDestroy(ErgmState *s);

#endif
