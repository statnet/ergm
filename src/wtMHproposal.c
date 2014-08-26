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
	     WtNetwork *nwp){

  char *fn, *sn;
  int i;
  for (i = 0; MHproposaltype[i] != ' ' && MHproposaltype[i] != 0; i++);
  MHproposaltype[i] = 0;
  /* Extract the required string information from the relevant sources */
  if((fn=(char *)malloc(sizeof(char)*(i+4)))==NULL){
    error("Error in MCMCSample: Can't allocate %d bytes for fn. Memory has not been deallocated, so restart R sometime soon.\n",
	  sizeof(char)*(i+4));
  }
  fn[0]='M';
  fn[1]='H';
  fn[2]='_';
  for(int j=0;j<i;j++)
    fn[j+3]=MHproposaltype[j];
  fn[i+3]='\0';
  /* fn is now the string 'MH_[name]', where [name] is MHproposaltype */
  for (i = 0; MHproposalpackage[i] != ' ' && MHproposalpackage[i] != 0; i++);
  MHproposalpackage[i] = 0;
  if((sn=(char *)malloc(sizeof(char)*(i+1)))==NULL){
    error("Error in ModelInitialize: Can't allocate %d bytes for sn. Memory has not been deallocated, so restart R sometime soon.\n",
	  sizeof(char)*(i+1));
  }
  sn=strncpy(sn,MHproposalpackage,i);
  sn[i]='\0';
  
  /* Search for the MH proposal function pointer */
  MHp->func=(void (*)(WtMHproposal*, WtNetwork*)) R_FindSymbol(fn,sn,NULL);
  if(MHp->func==NULL){
    error("Error in MH_* initialization: could not find function %s in "
	  "namespace for package %s."
	  "Memory has not been deallocated, so restart R sometime soon.\n",fn,sn);
  }

  MHp->inputs=inputs;
  
  MHp->discord=NULL;

  /*Clean up by freeing sn and fn*/
  free((void *)fn);
  free((void *)sn);

  MHp->ntoggles=0;
  (*(MHp->func))(MHp, nwp); /* Call MH proposal function to initialize */
  MHp->toggletail = (Vertex *)malloc(MHp->ntoggles * sizeof(Vertex));
  MHp->togglehead = (Vertex *)malloc(MHp->ntoggles * sizeof(Vertex));
  MHp->toggleweight = (double *)malloc(MHp->ntoggles * sizeof(double));
}

/*********************
 void WtMH_free

 A helper function to free memory allocated by WtMH_init.
*********************/
void WtMH_free(WtMHproposal *MHp){
  if(MHp->discord){
    for(WtNetwork **nwp=MHp->discord; *nwp!=NULL; nwp++){
      WtNetworkDestroy(*nwp);
    }
    free(MHp->discord);
  }
  free(MHp->toggletail);
  free(MHp->togglehead);
  free(MHp->toggleweight);
}

