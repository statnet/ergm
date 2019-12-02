#include "ergm_state.h"

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
                         Rboolean timings, int time, int *lasttoggle){

  ErgmState *s = Calloc(1, ErgmState);

  /* Form the network */
  s->nwp=NetworkInitialize(tails, heads, n_edges, 
                           n_nodes, directed_flag, bip, timings, time, lasttoggle);

  /* Initialize the model */
  s->m=ModelInitialize(funnames, sonames, &inputs, nterms, s->nwp, noinit_s);

  /* Initialize the M-H proposal */
  s->MHp=NULL;
  if(MHProposaltype)
    s->MHp = MHProposalInitialize(MHProposaltype, MHProposalpackage,
			       inputs,
			       s->nwp, attribs, maxout, maxin, minout, minin,
			       condAllDegExact, attriblength,
			       s->m->termarray->aux_storage);

  return s;
}

void ErgmStateDestroy(ErgmState *s){
  if(s->MHp) MHProposalDestroy(s->MHp, s->nwp);
  if(s->m) ModelDestroy(s->nwp, s->m);
  if(s->nwp) NetworkDestroy(s->nwp);
  Free(s);
}
