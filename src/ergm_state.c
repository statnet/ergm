#include "ergm_state.h"

ErgmState *ErgmStateInit(// Network settings
                         Vertex n_nodes, Rboolean directed_flag, Vertex bip,
                         // Model settings
                         SEXP mR, Rboolean noinit_s,
                         // Proposal settings
                         SEXP pR,
                         // Network state
                         Edge n_edges,
                         Vertex *tails, Vertex *heads,
                         Rboolean timings, int time, int *lasttoggle){

  ErgmState *s = Calloc(1, ErgmState);

  /* Form the network */
  s->nwp=NetworkInitialize(tails, heads, n_edges, 
                           n_nodes, directed_flag, bip, timings, time, lasttoggle);

  /* Initialize the model */
  s->m=NULL;
  if(s->nwp && mR) // Model also requires network.
    s->m = ModelInitialize(mR, s->nwp, noinit_s);

  /* Initialize the M-H proposal */
  s->MHp=NULL;
  if(s->m && pR) // Proposal also requires model's auxiliaries.
    s->MHp = MHProposalInitialize(pR, s->nwp, s->m->termarray->aux_storage);

  return s;
}

void ErgmStateDestroy(ErgmState *s){
  if(s->MHp) MHProposalDestroy(s->MHp, s->nwp);
  if(s->m) ModelDestroy(s->nwp, s->m);
  if(s->nwp) NetworkDestroy(s->nwp);
  Free(s);
}
