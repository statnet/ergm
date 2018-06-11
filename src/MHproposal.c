/*  File src/MHproposal.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2018 Statnet Commons
 */
#include "MHproposal.h"


/*********************
 void MH_init

 A helper function to process the MH_* related initialization.
*********************/
void MH_init(MHproposal *MHp, 
	     char *MHproposaltype, char *MHproposalpackage,
	     double *inputs,
	     int fVerbose,
	     Network *nwp,
	     int *attribs, int *maxout, int *maxin, 
	     int *minout, int *minin, int condAllDegExact, 
	     int attriblength){

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
  MHp->func=(void (*)(MHproposal*, Network*)) R_FindSymbol(fn,sn,NULL);
  if(MHp->func==NULL){
    error("Error in MH_* initialization: could not find function %s in "
	  "namespace for package %s."
	  "Memory has not been deallocated, so restart R sometime soon.\n",fn,sn);
  }

  MHp->inputs=inputs;

  MHp->bd=DegreeBoundInitialize(attribs, maxout, maxin, minout, minin,
			       condAllDegExact, attriblength, nwp);
  MHp->discord=NULL;

  /*Clean up by freeing sn and fn*/
  free((void *)fn);
  free((void *)sn);

  MHp->ntoggles=0;
  (*(MHp->func))(MHp, nwp); /* Call MH proposal function to initialize */
  MHp->toggletail = (Vertex *)malloc(MHp->ntoggles * sizeof(Vertex));
  MHp->togglehead = (Vertex *)malloc(MHp->ntoggles * sizeof(Vertex));
}

/*********************
 void MH_free

 A helper function to free memory allocated by MH_init.
*********************/
void MH_free(MHproposal *MHp){
  if(MHp->bd)DegreeBoundDestroy(MHp->bd);
  if(MHp->discord){
    for(Network **nwp=MHp->discord; *nwp!=NULL; nwp++){
      NetworkDestroy(*nwp);
    }
    free(MHp->discord);
  }
  free(MHp->toggletail);
  free(MHp->togglehead);
}

/***********************
 DegreeBound* DegreeBoundInitialize
************************/
DegreeBound* DegreeBoundInitialize(int *attribs, int *maxout, int *maxin,
	          	   int *minout, 
			   int *minin, int condAllDegExact,  int attriblength,
			   Network *nwp)
{
  int i,j;
  DegreeBound *bd;

  // This test no longer works, since the integer(0)->NULL no longer holds.
  // if(!(minout||minin||maxout||maxin||condAllDegExact)) return NULL; 
  
  if(!condAllDegExact && !attriblength) return NULL;
  

  bd = (DegreeBound *) malloc(sizeof(DegreeBound));

  bd->fBoundDegByAttr = 0;
  bd->attrcount = condAllDegExact ? 1 : attriblength / nwp->nnodes;
  bd->attribs = (int *) malloc(sizeof(int) * attriblength);
  bd->maxout  = (int *) malloc(sizeof(int) * attriblength);
  bd->maxin   = (int *) malloc(sizeof(int) * attriblength);
  bd->minout  = (int *) malloc(sizeof(int) * attriblength);
  bd->minin   = (int *) malloc(sizeof(int) * attriblength);
  
  /* bound by degree by attribute per node */
  if (bd->attrcount)
    {
      /* flag that we have data here */
      bd->fBoundDegByAttr = 1;
      
      if (!condAllDegExact)
	{
	  for (i=1; i <= nwp->nnodes; i++)
	    for (j=0; j < bd->attrcount; j++)
	      {
		bd->attribs[i-1 + j*nwp->nnodes] = 
		  attribs[(i - 1 + j*nwp->nnodes)];
		bd->maxout[i-1 + j*nwp->nnodes] =  
		  maxout[(i - 1 + j*nwp->nnodes)];
		bd->maxin[i-1 + j*nwp->nnodes] =  
		  maxin[(i - 1 + j*nwp->nnodes)];
		bd->minout[i-1 + j*nwp->nnodes] =  
		  minout[(i - 1 + j*nwp->nnodes)];
		bd->minin[i-1 + j*nwp->nnodes] =   
		  minin[(i - 1 + j*nwp->nnodes)];
	      }
	}
      else  /* condAllDegExact == TRUE */
	{
	  /* all ego columns get values of current in and out degrees;
	   max and min ego columns for (each of in and out) get same value; */
	  for (i=1;i<=nwp->nnodes;i++)
	    bd->maxout[i-1] = bd->minout[i-1] = nwp->outdegree[i];
	  
	  for (i=1;i<=nwp->nnodes;i++)
	    bd->maxin[i-1] = bd->minin[i-1] = nwp->indegree[i];
	}
      return bd;
    }
  else
    {
      return NULL;
    }
}


/*****************
  void DegreeBoundDestroy
******************/
void DegreeBoundDestroy(DegreeBound *bd)
{  
  free(bd->attribs); 
  free(bd->maxout); 
  free(bd->minout); 
  free(bd->maxin); 
  free(bd->minin); 
  free(bd);
}


/********************
 int CheckTogglesValid
********************/
int CheckTogglesValid(MHproposal *MHp, Network *nwp) {
  int fvalid;
  int i;
  DegreeBound *bd=MHp->bd;

  if(!bd) return 1;

  /* *** don't forget when getting attributes that tail-> head */
  int *tailattr = (int *) malloc(sizeof(int) * bd->attrcount);
  int *headattr = (int *) malloc(sizeof(int) * bd->attrcount);
  
  fvalid = 1;
  
  /* Make proposed toggles */
  for (i=0; i<MHp->ntoggles; i++)
    ToggleEdge(MHp->toggletail[i], MHp->togglehead[i], nwp);

  /*  Rprintf("fvalid %d bd->fBoundDegByAttr %d\n", fvalid, bd->fBoundDegByAttr); */

  /* if we're bounding degrees by attribute */
  if (bd->fBoundDegByAttr && fvalid) {
    Edge e;
    Vertex v;
    int k; 
    if (nwp->directed_flag) {
      /* for each tail and head pair */
      for (i = 0; i < MHp->ntoggles && fvalid; i++) {
        /* work through each attribute for each toggle */
        for (k=0; k < bd->attrcount; k++){
	        tailattr[k] = headattr[k] = 0;
        }
        /* calculate tail outdegree totals for each attribute
        for each outedge of the tail 	      */
	      
        for(e = EdgetreeMinimum(nwp->outedges, MHp->toggletail[i]);
        (v = nwp->outedges[e].value) != 0;
        e = EdgetreeSuccessor(nwp->outedges, e)) {
          for (k=0; k < bd->attrcount; k++)
            if (bd->attribs[v-1 + k*nwp->nnodes]) tailattr[k]++;
        }
	      
        /* calculate head indegree totals for each attribute
        for each inedge of the head */
	      
        for(e = EdgetreeMinimum(nwp->inedges, MHp->togglehead[i]);
        (v = nwp->inedges[e].value) != 0;
        e = EdgetreeSuccessor(nwp->inedges, e)) {
          for (k=0; k < bd->attrcount; k++)
            if (bd->attribs[v-1 + k*nwp->nnodes]) headattr[k]++;
        }

        /* for each attribute */

        for (k=0; k < bd->attrcount && fvalid; k++){
          fvalid=!((tailattr[k]>bd->maxout[MHp->toggletail[i]-1+k*nwp->nnodes])||
          (tailattr[k] < bd->minout[MHp->toggletail[i]-1+k*nwp->nnodes]) || 
          (headattr[k] >  bd->maxin[MHp->togglehead[i]-1+k*nwp->nnodes]) ||
          (headattr[k] <  bd->minin[MHp->togglehead[i]-1+k*nwp->nnodes]) );
        }
      }
    }
    else { /* ! nwp->directed_flag  */
      /* for each tail and head pair, (in that order: (tail, head)) */
      for (i = 0; i < MHp->ntoggles && fvalid; i++) {
        for (k=0; k < bd->attrcount; k++){
	        tailattr[k] = headattr[k] = 0;
	      }
	      
	      /* calculate tail totals for each attribute
        for each outedge and inedge of the tail  */
	      
	      for(e = EdgetreeMinimum(nwp->outedges, MHp->toggletail[i]);
        (v = nwp->outedges[e].value) != 0;
        e = EdgetreeSuccessor(nwp->outedges, e)) {
          for (k=0; k < bd->attrcount; k++)
            if (bd->attribs[v-1 + k*nwp->nnodes])
              tailattr[k]++;
        }
	      for(e = EdgetreeMinimum(nwp->inedges, MHp->toggletail[i]);
        (v = nwp->inedges[e].value) != 0;
        e = EdgetreeSuccessor(nwp->inedges, e)) {
          for (k=0; k < bd->attrcount; k++)
            if (bd->attribs[v-1 + k*nwp->nnodes])
              tailattr[k]++;
        }
	      
	      /* calculate head totals for each attribute
        for each outedge and inedge of the head */
	      
	      for(e = EdgetreeMinimum(nwp->outedges, MHp->togglehead[i]);
        (v = nwp->outedges[e].value) != 0;
        e = EdgetreeSuccessor(nwp->outedges, e)) {
          for (k=0; k < bd->attrcount; k++)
            if (bd->attribs[v-1 + k*nwp->nnodes])
              headattr[k]++;
        }
	      for(e = EdgetreeMinimum(nwp->inedges, MHp->togglehead[i]);
        (v = nwp->inedges[e].value) != 0;
        e = EdgetreeSuccessor(nwp->inedges, e)) {
          for (k=0; k < bd->attrcount; k++)
            if (bd->attribs[v-1 + k*nwp->nnodes])
              headattr[k]++;
        }

	      /* for each attribute
        check tails' and heads' outmax and outmin */
	      for (k=0; k < bd->attrcount && fvalid; k++){
          fvalid=!((tailattr[k]>bd->maxout[MHp->toggletail[i]-1+k*nwp->nnodes])|| 
          (tailattr[k] < bd->minout[MHp->toggletail[i]-1+k*nwp->nnodes]) || 
          (headattr[k] > bd->maxout[MHp->togglehead[i]-1+k*nwp->nnodes]) ||
          (headattr[k] < bd->minout[MHp->togglehead[i]-1+k*nwp->nnodes]) );
	      }
	    }
    }
  }
  
  free(tailattr);
  free(headattr);
  
  /* Undo proposed toggles (of edges(tail, head)) */
  for (i=0; i<MHp->ntoggles; i++)
    ToggleEdge(MHp->toggletail[i], MHp->togglehead[i], nwp);

  return fvalid;
}

int CheckConstrainedTogglesValid(MHproposal *MHp, Network *nwp)
{
  int fvalid = 1;
  int i;
  DegreeBound *bd=MHp->bd;

  if(!bd) return 1;

  /* Make proposed toggles */
  for (i=0; i<MHp->ntoggles; i++)
    ToggleEdge(MHp->toggletail[i], MHp->togglehead[i], nwp);

  /* if we're bounding degrees by attribute */
  if (bd->fBoundDegByAttr && fvalid)
  {
    Edge e;
    Vertex v;
    int k;
    int *tailattr = (int *) malloc(sizeof(int) * bd->attrcount);
    int *headattr = (int *) malloc(sizeof(int) * bd->attrcount);
    
    if (nwp->directed_flag)
    {
      /* for each tail and head pair - yes (tail, head), not (head,tail) */
      for (i = 0; i < MHp->ntoggles && fvalid; i++) {
        for (k=0; k < bd->attrcount; k++){
	        tailattr[k] = headattr[k] = 0;
	      }
	      /* calculate tail outdegree totals for each attribute
        for each outedge of the tail 	      */
	      
	      for(e = EdgetreeMinimum(nwp->outedges, MHp->toggletail[i]);
        (v = nwp->outedges[e].value) != 0;
        e = EdgetreeSuccessor(nwp->outedges, e))
        {
          for (k=0; k < bd->attrcount; k++)
            if (bd->attribs[v-1 + k*nwp->nnodes])
              tailattr[k]++;
        }
	      
	      /* calculate head indegree totals for each attribute
        for each inedge of the head */
	      
	      for(e = EdgetreeMinimum(nwp->inedges, MHp->togglehead[i]);
        (v = nwp->inedges[e].value) != 0;
        e = EdgetreeSuccessor(nwp->inedges, e))
        {
          for (k=0; k < bd->attrcount; k++)
            if (bd->attribs[v-1 + k*nwp->nnodes])
              headattr[k]++;
        }

	      /* for each attribute */
	      for (k=0; k < bd->attrcount && fvalid; k++){
          fvalid=!((tailattr[k]>bd->maxout[MHp->toggletail[i]-1+k*nwp->nnodes])||
          (tailattr[k] < bd->minout[MHp->toggletail[i]-1+k*nwp->nnodes]) || 
          (headattr[k] >  bd->maxin[MHp->togglehead[i]-1+k*nwp->nnodes]) ||
          (headattr[k] <  bd->minin[MHp->togglehead[i]-1+k*nwp->nnodes])) ;
	      }
	    }
    }
    else /* ! nwp->directed_flag */
    {
      /* for each tail and head pair */
      for (i = 0; i < MHp->ntoggles && fvalid; i++)
	    {
        for (k=0; k < bd->attrcount; k++){
	        tailattr[k] = headattr[k] = 0;
	      }
	      
	      /* calculate tail totals for each attribute
        for each outedge and inedge of the tail  */
	      
	      for(e = EdgetreeMinimum(nwp->outedges, MHp->toggletail[i]);
        (v = nwp->outedges[e].value) != 0;
        e = EdgetreeSuccessor(nwp->outedges, e))
        {
          for (k=0; k < bd->attrcount; k++)
            if (bd->attribs[v-1 + k*nwp->nnodes])
              tailattr[k]++;
        }
	      for(e = EdgetreeMinimum(nwp->inedges, MHp->toggletail[i]);
        (v = nwp->inedges[e].value) != 0;
        e = EdgetreeSuccessor(nwp->inedges, e))
        {
          for (k=0; k < bd->attrcount; k++)
            if (bd->attribs[v-1 + k*nwp->nnodes])
              tailattr[k]++;
        }
	      
	      /* calculate head totals for each attribute
        for each outedge and inedge of the head */
	      for(e = EdgetreeMinimum(nwp->outedges, MHp->togglehead[i]);
        (v = nwp->outedges[e].value) != 0;
        e = EdgetreeSuccessor(nwp->outedges, e))
        {
          for (k=0; k < bd->attrcount; k++)
            if (bd->attribs[v-1 + k*nwp->nnodes])
              headattr[k]++;
        }
	      for(e = EdgetreeMinimum(nwp->inedges, MHp->togglehead[i]);
        (v = nwp->inedges[e].value) != 0;
        e = EdgetreeSuccessor(nwp->inedges, e))
        {
          for (k=0; k < bd->attrcount; k++)
            if (bd->attribs[v-1 + k*nwp->nnodes])
              headattr[k]++;
        }
        
	      /* for each attribute
        check tails' and heads' outmax and outmin */
	      for (k=0; k < bd->attrcount && fvalid; k++)
          fvalid=!(tailattr[k]>bd->maxout[MHp->toggletail[i]-1+k*nwp->nnodes])||
        (tailattr[k] < bd->minout[MHp->toggletail[i]-1+k*nwp->nnodes]) || 
        (headattr[k] > bd->maxout[MHp->togglehead[i]-1+k*nwp->nnodes]) ||
        (headattr[k] < bd->minout[MHp->togglehead[i]-1+k*nwp->nnodes]) ;
	    }
    }
    free(tailattr);
    free(headattr);
  }
  /* Make proposed toggles (of edges (tail, head), not (head, tail) */
  for (i=0; i<MHp->ntoggles; i++)
    ToggleEdge(MHp->toggletail[i], MHp->togglehead[i], nwp);
  
  return fvalid;
}
