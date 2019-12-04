#include "ergm_wtstate.h"

WtErgmState *WtErgmStateInit(// Network settings
			     Vertex n_nodes, Rboolean directed_flag, Vertex bip,
                             // Model settings
			     SEXP mR, Rboolean noinit_s,
                             // Proposal settings
			     SEXP pR,
                             // Network state
                             Edge n_edges,
			     Vertex *tails, Vertex *heads, double *weights,
                             Rboolean timings, int time, int *lasttoggle){
  WtErgmState *s = Calloc(1, WtErgmState);
  
  /* Form the network */
  s->nwp=WtNetworkInitialize(tails, heads, weights, n_edges, 
                             n_nodes, directed_flag, bip, timings, time, lasttoggle);

  /* Initialize the model */
  s->m=NULL;
  if(s->nwp && mR) // Model also requires network.
    s->m = WtModelInitialize(mR, s->nwp, noinit_s);

  /* Initialize the M-H proposal */
  s->MHp=NULL;
  if(s->m && pR) // Proposal also requires model's auxiliaries.
    s->MHp = WtMHProposalInitialize(pR, s->nwp, s->m->termarray->aux_storage);
  
  return s;
}

void WtErgmStateDestroy(WtErgmState *s){
  if(s->MHp) WtMHProposalDestroy(s->MHp, s->nwp);
  if(s->m) WtModelDestroy(s->nwp, s->m);
  if(s->nwp) WtNetworkDestroy(s->nwp);
  Free(s);
}
