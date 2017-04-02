/*  File src/wtMHproposal.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2013 Statnet Commons
 */
#include "wtMHproposal.h"


/*********************
 void WtMH_init

 A helper function to process the MH_* related initialization.
*********************/
void WtMH_init(WtMHproposal *MHp, 
	     char *MHproposaltype, char *MHproposalpackage, 
	       double *inputs,
	     int fVerbose,
	     WtNetwork *nwp,
	     void **aux_storage){

  char *fn, *sn;
  int i;

  MHp->i_func=MHp->p_func=MHp->f_func=NULL;
  MHp->u_func=NULL;
  MHp->storage=NULL;
  
  for (i = 0; MHproposaltype[i] != ' ' && MHproposaltype[i] != 0; i++);
  MHproposaltype[i] = 0;
  /* Extract the required string information from the relevant sources */
  fn=Calloc(i+4, char);
  fn[0]='M';
  fn[1]='H';
  fn[2]='_';
  for(int j=0;j<i;j++)
    fn[j+3]=MHproposaltype[j];
  fn[i+3]='\0';
  /* fn is now the string 'MH_[name]', where [name] is MHproposaltype */
  for (i = 0; MHproposalpackage[i] != ' ' && MHproposalpackage[i] != 0; i++);
  MHproposalpackage[i] = 0;
  sn=Calloc(i+1, char);
  sn=strncpy(sn,MHproposalpackage,i);
  sn[i]='\0';
  
  /* Search for the MH proposal function pointer */
  // Old-style name:
  MHp->p_func=(void (*)(WtMHproposal*, WtNetwork*)) R_FindSymbol(fn,sn,NULL);
  if(MHp->p_func==NULL){
    // New-style name:
    fn[1] = 'p';
    MHp->p_func=(void (*)(WtMHproposal*, WtNetwork*)) R_FindSymbol(fn,sn,NULL);
    if(MHp->p_func==NULL){    
      error("Error in the proposal initialization: could not find function %s in "
	    "namespace for package %s."
	    "Memory has not been deallocated, so restart R sometime soon.\n",fn,sn);
    }
  }

  // Optional functions
  fn[1] = 'i';
  MHp->i_func=(void (*)(WtMHproposal*, WtNetwork*)) R_FindSymbol(fn,sn,NULL);
  fn[1] = 'u';
  MHp->u_func=(void (*)(Vertex tail, Vertex head, double weight, WtMHproposal*, WtNetwork*)) R_FindSymbol(fn,sn,NULL);
  fn[1] = 'f';
  MHp->f_func=(void (*)(WtMHproposal*, WtNetwork*)) R_FindSymbol(fn,sn,NULL);
    
  MHp->inputs=inputs;

  /*Clean up by freeing sn and fn*/
  Free(fn);
  Free(sn);

  MHp->aux_storage = aux_storage;

  MHp->ntoggles=0;
  if(MHp->i_func){
    // New-style initialization
    (*(MHp->i_func))(MHp, nwp);
  }else{
    // Old-style initialization
    (*(MHp->p_func))(MHp, nwp); /* Call MH proposal function to initialize */
  }
  
  MHp->toggletail = (Vertex *)Calloc(MHp->ntoggles, Vertex);
  MHp->togglehead = (Vertex *)Calloc(MHp->ntoggles, Vertex);
  MHp->toggleweight = (double *)Calloc(MHp->ntoggles, double);
}

/*********************
 void WtMH_free

 A helper function to free memory allocated by WtMH_init.
*********************/
void WtMH_free(WtMHproposal *MHp, WtNetwork *nwp){
  if(MHp->f_func) (*(MHp->f_func))(MHp, nwp);
  if(MHp->storage){
    Free(MHp->storage);
    MHp->storage=NULL;
  }
  MHp->aux_storage=NULL;
  Free(MHp->toggletail);
  Free(MHp->togglehead);
  Free(MHp->toggleweight);
}

