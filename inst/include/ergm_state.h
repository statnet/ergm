#ifndef _ERGM_STATE_H_
#define _ERGM_STATE_H_

#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_MHproposal.h"
#include "ergm_model.h"

#define ARGS_NWSETTINGS SEXP dn, SEXP dflag, SEXP bipartite
#define ARGS_MODEL SEXP nterms, SEXP funnames, SEXP sonames
#define ARGS_MHPROPOSAL SEXP MHProposaltype, SEXP MHProposalpackage,    \
    SEXP attribs, SEXP maxout, SEXP maxin, SEXP minout,                 \
    SEXP minin, SEXP condAllDegExact, SEXP attriblength
#define ARGS_INPUTS SEXP inputs
#define ARGS_NWSTATE SEXP nedges, SEXP tails, SEXP heads
#define ARGS_LASTTOGGLE SEXP time, SEXP lasttoggle

#define YES_NWSETTINGS asInteger(dn), asInteger(dflag), asInteger(bipartite)
#define YES_MODEL asInteger(nterms), FIRSTCHAR(funnames), FIRSTCHAR(sonames), FALSE
#define YES_MODEL_NOINIT_S asInteger(nterms), FIRSTCHAR(funnames), FIRSTCHAR(sonames), TRUE
#define YES_MHPROPOSAL FIRSTCHAR(MHProposaltype), FIRSTCHAR(MHProposalpackage), \
    INTEGER(attribs), INTEGER(maxout), INTEGER(maxin),                  \
    INTEGER(minout), INTEGER(minin), asInteger(condAllDegExact), asInteger(attriblength)
#define YES_INPUTS REAL(inputs)
#define YES_NWSTATE asInteger(nedges), (Vertex*) INTEGER(tails), (Vertex*) INTEGER(heads)
#define YES_LASTTOGGLE length(time)>0, asInteger(time), INTEGER(lasttoggle)

#define NO_MODEL 0, NULL, NULL, FALSE
#define NO_MHPROPOSAL  NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0
#define NO_INPUTS NULL
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
