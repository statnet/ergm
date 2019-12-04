#ifndef _ERGM_WTSTATE_H_
#define _ERGM_WTSTATE_H_

#include "ergm_wtedgetree.h"
#include "ergm_wtchangestat.h"
#include "ergm_wtMHproposal.h"
#include "ergm_wtmodel.h"

#define ARGS_WTNWSETTINGS SEXP dn, SEXP dflag, SEXP bipartite
#define ARGS_WTMODEL SEXP mR
#define ARGS_WTMHPROPOSAL SEXP pR
#define ARGS_WTNWSTATE SEXP nedges, SEXP tails, SEXP heads, SEXP weights
#define ARGS_WTLASTTOGGLE SEXP time, SEXP lasttoggle

#define YES_WTNWSETTINGS asInteger(dn), asInteger(dflag), asInteger(bipartite)
#define YES_WTMODEL mR, FALSE
#define YES_WTMODEL_NOINIT_S mR, TRUE
#define YES_WTMHPROPOSAL pR
#define YES_WTNWSTATE asInteger(nedges), (Vertex*) INTEGER(tails), (Vertex*) INTEGER(heads), REAL(weights)
#define YES_WTLASTTOGGLE length(time)>0, asInteger(time), INTEGER(lasttoggle)

#define NO_WTMODEL NULL, FALSE
#define NO_WTMHPROPOSAL  NULL
#define NO_WTLASTTOGGLE FALSE, 0, NULL
#define NO_WTNWSTATE 0, NULL, NULL, NULL

typedef struct{
  WtNetwork *nwp;
  WtModel *m;
  WtMHProposal *MHp;
} WtErgmState;

WtErgmState *WtErgmStateInit(// Network settings
			     Vertex n_nodes, Rboolean directed_flag, Vertex bip,
                             // Model settings
			     SEXP mR, Rboolean noinit_s,
                             // Proposal settings
                             SEXP pR,
                             // Network state
                             Edge n_edges,
			     Vertex *tails, Vertex *heads, double *weights,
                             Rboolean timings, int time, int *lasttoggle);
void WtErgmStateDestroy(WtErgmState *s);

#endif
