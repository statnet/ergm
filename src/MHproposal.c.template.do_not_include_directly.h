/*  File src/MHproposal.c.template.do_not_include_directly.h in package ergm,
 *  part of the Statnet suite of packages for network analysis,
 *  https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */
#include "ergm_Rutil.h"

void ETYPE(On,NetworkEdgeChangeMUWrap)(Vertex tail, Vertex head, IFEWT(EWTTYPE weight,) void *MHp, ETYPE(Network) *nwp, EWTTYPE edgestate){
  ((ETYPE(MHProposal) *) MHp)->u_func(tail, head, IFEWT(weight,) MHp, nwp, edgestate);
}

/*********************
 void ETYPE(MHProposalInitialize)

 A helper function to process the MH_* related initialization.
*********************/
ETYPE(MHProposal) *ETYPE(MHProposalInitialize)(SEXP pR, ETYPE(Network) *nwp, void **aux_storage){
  ETYPE(MHProposal) *MHp = R_Calloc(1, ETYPE(MHProposal));
  MHp->R = pR;

  MHp->i_func=MHp->p_func=MHp->f_func=NULL;
  MHp->u_func=NULL;
  MHp->storage=NULL;
  
  /* Extract the required string information from the relevant sources */
  const char *fname = FIRSTCHAR(getListElement(pR, "name")),
    *sn = FIRSTCHAR(getListElement(pR, "pkgname"));

  char *fn = R_Calloc(strlen(fname)+4, char);
  fn[0]='M';
  fn[1]='H';
  fn[2]='_';
  strcpy(fn+3, fname);
  
  /* Search for the MH proposal function pointer */
  // Old-style name:
  MHp->p_func=(void (*)(ETYPE(MHProposal)*, ETYPE(Network)*)) R_FindSymbol(fn,sn,NULL);
  if(MHp->p_func==NULL){
    // New-style name:
    fn[1] = 'p';
    MHp->p_func=(void (*)(ETYPE(MHProposal)*, ETYPE(Network)*)) R_FindSymbol(fn,sn,NULL);
    if(MHp->p_func==NULL){    
      error("Error in the proposal initialization: could not find function %s in "
	  "namespace for package %s."
	  "Memory has not been deallocated, so restart R sometime soon.\n",fn,sn);
    }
  }

  // Optional functions
  fn[1] = 'i';
  MHp->i_func=(void (*)(ETYPE(MHProposal)*, ETYPE(Network)*)) R_FindSymbol(fn,sn,NULL);
  fn[1] = 'u';
  MHp->u_func=(void (*)(Vertex, Vertex, IFEWT(EWTTYPE,) ETYPE(MHProposal)*, ETYPE(Network)*, EWTTYPE)) R_FindSymbol(fn,sn,NULL);
  fn[1] = 'f';
  MHp->f_func=(void (*)(ETYPE(MHProposal)*, ETYPE(Network)*)) R_FindSymbol(fn,sn,NULL);
  fn[1] = 'x';
  MHp->x_func=(void (*)(unsigned int, void *, ETYPE(MHProposal)*, ETYPE(Network)*)) R_FindSymbol(fn,sn,NULL);

  SEXP tmp = getListElement(pR, "inputs");
  MHp->ninputs = length(tmp);
  MHp->inputs = MHp->ninputs ? REAL(tmp) : NULL;
  tmp = getListElement(pR, "iinputs");
  MHp->niinputs = length(tmp);
  MHp->iinputs = MHp->niinputs ? INTEGER(tmp) : NULL;

  /*Clean up by freeing sn and fn*/
  R_Free(fn);

  MHp->aux_storage = aux_storage;
  SEXP aux_slotsR = getListElement(pR,"aux.slots");
  if((MHp->n_aux = length(aux_slotsR))){
    MHp->aux_slots = (unsigned int *) INTEGER(aux_slotsR);
  }else MHp->aux_slots = NULL;
  
  MHp->ntoggles=0;
  if(MHp->i_func){
    // New-style initialization
    (*(MHp->i_func))(MHp, nwp);
  }else{
    // Old-style initialization
    (*(MHp->p_func))(MHp, nwp); /* Call MH proposal function to initialize */
  }
  
  if(MHp->ntoggles==MH_FAILED){
    REprintf("MH proposal function's initial network configuration is one from which no toggle(s) can be proposed.\n");
    MHp->toggletail = MHp->togglehead = NULL; // To be safe.
    IFEWT(MHp->toggleweight = NULL);
    MHp->u_func = NULL; // Important, since the callback was never installed.
    ETYPE(MHProposalDestroy)(MHp, nwp);
    return NULL;
  }
  MHp->toggletail = (Vertex *)R_Calloc(MHp->ntoggles, Vertex);
  MHp->togglehead = (Vertex *)R_Calloc(MHp->ntoggles, Vertex);
  IFEWT(MHp->toggleweight = (EWTTYPE *)R_Calloc(MHp->ntoggles, EWTTYPE));

  if(MHp->u_func){
    ETYPE(AddOn,NetworkEdgeChange)(nwp, ETYPE(On,NetworkEdgeChangeMUWrap), MHp, 0); // Need to insert at the start.
  }

  return MHp;
}

/*********************
 void ETYPE(MHProposalDestroy)

 A helper function to free memory allocated by ETYPE(MHProposalInitialize).
*********************/
void ETYPE(MHProposalDestroy)(ETYPE(MHProposal) *MHp, ETYPE(Network) *nwp){
  if(!MHp) return;
  if(MHp->u_func) ETYPE(DeleteOn,NetworkEdgeChange)(nwp, ETYPE(On,NetworkEdgeChangeMUWrap), MHp);
  if(MHp->f_func) (*(MHp->f_func))(MHp, nwp);
  if(MHp->storage){
    R_Free(MHp->storage);
    MHp->storage=NULL;
  }
  MHp->aux_storage=NULL;
  R_Free(MHp->toggletail);
  R_Free(MHp->togglehead);
  IFEWT(R_Free(MHp->toggleweight));

  R_Free(MHp);
}

