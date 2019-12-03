#ifndef _ERGM_WTSTATE_H_
#define _ERGM_WTSTATE_H_

#include "ergm_wtedgetree.h"
#include "ergm_wtchangestat.h"
#include "ergm_wtMHproposal.h"
#include "ergm_wtmodel.h"

#define ARGS_WTNWSETTINGS SEXP dn, SEXP dflag, SEXP bipartite
#define ARGS_WTMODEL SEXP nterms, SEXP funnames, SEXP sonames
#define ARGS_WTMHPROPOSAL SEXP MHProposaltype, SEXP MHProposalpackage
#define ARGS_WTINPUTS SEXP inputs
#define ARGS_WTNWSTATE SEXP nedges, SEXP tails, SEXP heads, SEXP weights
#define ARGS_WTLASTTOGGLE SEXP time, SEXP lasttoggle

#define YES_WTNWSETTINGS asInteger(dn), asInteger(dflag), asInteger(bipartite)
#define YES_WTMODEL asInteger(nterms), FIRSTCHAR(funnames), FIRSTCHAR(sonames), FALSE
#define YES_WTMODEL_NOINIT_S asInteger(nterms), FIRSTCHAR(funnames), FIRSTCHAR(sonames), TRUE
#define YES_WTMHPROPOSAL FIRSTCHAR(MHProposaltype), FIRSTCHAR(MHProposalpackage)
#define YES_WTINPUTS REAL(inputs)
#define YES_WTNWSTATE asInteger(nedges), (Vertex*) INTEGER(tails), (Vertex*) INTEGER(heads), REAL(weights)
#define YES_WTLASTTOGGLE length(time)>0, asInteger(time), INTEGER(lasttoggle)

#define NO_WTMODEL 0, NULL, NULL, FALSE
#define NO_WTMHPROPOSAL  NULL, NULL
#define NO_WTINPUTS NULL
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
			     int nterms, const char *funnames, const char *sonames, Rboolean noinit_s,
                             // Proposal settings
			     const char *MHProposaltype, const char *MHProposalpackage,
			     double *inputs,
                             // Network state
                             Edge n_edges,
			     Vertex *tails, Vertex *heads, double *weights,
                             Rboolean timings, int time, int *lasttoggle);
void WtErgmStateDestroy(WtErgmState *s);

#endif
