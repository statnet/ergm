#include "wtMHproposal.h"


/*********************
 void WtMH_init

 A helper function to process the MH_* related initialization.
*********************/
void WtMH_init(WtMHproposal *MH, 
	     char *MHproposaltype, char *MHproposalpackage, 
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
  MH->func=(void (*)(WtMHproposal*, WtNetwork*)) R_FindSymbol(fn,sn,NULL);
  if(MH->func==NULL){
    error("Error in MH_* initialization: could not find function %s in "
	  "namespace for package %s."
	  "Memory has not been deallocated, so restart R sometime soon.\n",fn,sn);
  }      
  
  /*Clean up by freeing sn and fn*/
  free((void *)fn);
  free((void *)sn);

  MH->ntoggles=0;
  (*(MH->func))(MH, nwp); /* Call MH proposal function to initialize */
  MH->togglehead = (Vertex *)malloc(MH->ntoggles * sizeof(Vertex));
  MH->toggletail = (Vertex *)malloc(MH->ntoggles * sizeof(Vertex));
  MH->toggleweight = (double *)malloc(MH->ntoggles * sizeof(double));
  MH->oldweight = (double *)malloc(MH->ntoggles * sizeof(double));
}

/*********************
 void WtMH_free

 A helper function to free memory allocated by WtMH_init.
*********************/
void WtMH_free(WtMHproposal *MH){
  free(MH->togglehead);
  free(MH->toggletail);
  free(MH->toggleweight);
  free(MH->oldweight);
}

