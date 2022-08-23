/*  File src/wtMHproposal.c in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2022 Statnet Commons
 */
#include "ergm_wtMHproposal.h"
#include "ergm_Rutil.h"

/*********************
 void WtMHProposalInitialize

 A helper function to process the MH_* related initialization.
*********************/
WtMHProposal *WtMHProposalInitialize(SEXP pR, WtNetwork *nwp, void **aux_storage){
  WtMHProposal *MHp = Calloc(1, WtMHProposal);
  MHp->R = pR;

  MHp->i_func=MHp->p_func=MHp->f_func=NULL;
  MHp->u_func=NULL;
  MHp->storage=NULL;
  
  /* Extract the required string information from the relevant sources */
  const char *fname = FIRSTCHAR(getListElement(pR, "name")),
    *sn = FIRSTCHAR(getListElement(pR, "pkgname"));
  char *fn = Calloc(strlen(fname)+4, char);
  fn[0]='M';
  fn[1]='H';
  fn[2]='_';
  strcpy(fn+3, fname);
  
  /* Search for the MH proposal function pointer */
  // Old-style name:
  MHp->p_func=(void (*)(WtMHProposal*, WtNetwork*)) R_FindSymbol(fn,sn,NULL);
  if(MHp->p_func==NULL){
    // New-style name:
    fn[1] = 'p';
    MHp->p_func=(void (*)(WtMHProposal*, WtNetwork*)) R_FindSymbol(fn,sn,NULL);
    if(MHp->p_func==NULL){    
      error("Error in the proposal initialization: could not find function %s in "
	  "namespace for package %s."
	  "Memory has not been deallocated, so restart R sometime soon.\n",fn,sn);
    }
  }

  // Optional functions
  fn[1] = 'i';
  MHp->i_func=(void (*)(WtMHProposal*, WtNetwork*)) R_FindSymbol(fn,sn,NULL);
  fn[1] = 'u';
  MHp->u_func=(void (*)(Vertex tail, Vertex head, double weight, WtMHProposal*, WtNetwork*)) R_FindSymbol(fn,sn,NULL);
  fn[1] = 'f';
  MHp->f_func=(void (*)(WtMHProposal*, WtNetwork*)) R_FindSymbol(fn,sn,NULL);
  fn[1] = 'x';
  MHp->x_func=(void (*)(unsigned int type, void *data, WtMHProposal*, WtNetwork*)) R_FindSymbol(fn,sn,NULL);

  SEXP tmp = getListElement(pR, "inputs");
  MHp->inputs=length(tmp) ? REAL(tmp) : NULL;
  tmp = getListElement(pR, "iinputs");
  MHp->iinputs=length(tmp) ? INTEGER(tmp) : NULL;

  /*Clean up by freeing sn and fn*/
  Free(fn);

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
    MHp->toggleweight = NULL;
    MHp->u_func = NULL; // Important, since the callback was never installed.
    WtMHProposalDestroy(MHp, nwp);
    return NULL;
  }
  MHp->toggletail = (Vertex *)Calloc(MHp->ntoggles, Vertex);
  MHp->togglehead = (Vertex *)Calloc(MHp->ntoggles, Vertex);
  MHp->toggleweight = (double *)Calloc(MHp->ntoggles, double);

  if(MHp->u_func){
    AddOnWtNetworkEdgeChange(nwp, (OnWtNetworkEdgeChange) MHp->u_func, MHp, 0); // Need to insert at the start.
  }

  return MHp;
}

/*********************
 void WtMHProposalDestroy

 A helper function to free memory allocated by WtMHProposalInitialize.
*********************/
void WtMHProposalDestroy(WtMHProposal *MHp, WtNetwork *nwp){
  if(!MHp) return;
  if(MHp->u_func) DeleteOnWtNetworkEdgeChange(nwp, (OnWtNetworkEdgeChange) MHp->u_func, MHp);
  if(MHp->f_func) (*(MHp->f_func))(MHp, nwp);
  if(MHp->storage){
    Free(MHp->storage);
    MHp->storage=NULL;
  }
  MHp->aux_storage=NULL;
  Free(MHp->toggletail);
  Free(MHp->togglehead);
  Free(MHp->toggleweight);

  Free(MHp);
}

