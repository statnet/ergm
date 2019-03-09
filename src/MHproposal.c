/*  File src/MHproposal.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2019 Statnet Commons
 */
#include "ergm_MHproposal.h"


/*********************
 void MHProposalInitialize

 A helper function to process the MH_* related initialization.
*********************/
MHProposal *MHProposalInitialize(
	     char *MHProposaltype, char *MHProposalpackage,
	     double *inputs,
	     int fVerbose,
	     Network *nwp,
	     int *attribs, int *maxout, int *maxin, 
	     int *minout, int *minin, int condAllDegExact, 
	     int attriblength){
  MHProposal *MHp = Calloc(1, MHProposal);
  
  char *fn, *sn;
  int i;
  for (i = 0; MHProposaltype[i] != ' ' && MHProposaltype[i] != 0; i++);
  MHProposaltype[i] = 0;
  /* Extract the required string information from the relevant sources */
  fn = Calloc(i+4, char);
  fn[0]='M';
  fn[1]='H';
  fn[2]='_';
  for(int j=0;j<i;j++)
    fn[j+3]=MHProposaltype[j];
  fn[i+3]='\0';
  /* fn is now the string 'MH_[name]', where [name] is MHProposaltype */
  for (i = 0; MHProposalpackage[i] != ' ' && MHProposalpackage[i] != 0; i++);
  MHProposalpackage[i] = 0;
  sn = Calloc(i+1, char);
  sn=strncpy(sn,MHProposalpackage,i);
  sn[i]='\0';
  
  /* Search for the MH proposal function pointer */
  MHp->func=(void (*)(MHProposal*, Network*)) R_FindSymbol(fn,sn,NULL);
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
  Free(fn);
  Free(sn);

  MHp->ntoggles=0;
  (*(MHp->func))(MHp, nwp); /* Call MH proposal function to initialize */
  MHp->toggletail = (Vertex *)Calloc(MHp->ntoggles, Vertex);
  MHp->togglehead = (Vertex *)Calloc(MHp->ntoggles, Vertex);

  return MHp;
}

/*********************
 void MHProposalDestroy

 A helper function to free memory allocated by MHProposalInitialize.
*********************/
void MHProposalDestroy(MHProposal *MHp){
  if(MHp->bd)DegreeBoundDestroy(MHp->bd);
  if(MHp->discord){
    for(Network **nwp=MHp->discord; *nwp!=NULL; nwp++){
      NetworkDestroy(*nwp);
    }
    Free(MHp->discord);
  }
  Free(MHp->toggletail);
  Free(MHp->togglehead);

  Free(MHp);
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
  

  bd = (DegreeBound *) Calloc(1, DegreeBound);

  bd->fBoundDegByAttr = 0;
  bd->attrcount = condAllDegExact ? 1 : attriblength / nwp->nnodes;
  bd->attribs = (int *) Calloc(attriblength, int);
  bd->maxout  = (int *) Calloc(attriblength, int);
  bd->maxin   = (int *) Calloc(attriblength, int);
  bd->minout  = (int *) Calloc(attriblength, int);
  bd->minin   = (int *) Calloc(attriblength, int);
  
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
  Free(bd->attribs); 
  Free(bd->maxout); 
  Free(bd->minout); 
  Free(bd->maxin); 
  Free(bd->minin); 
  Free(bd);
}


/********************
 int CheckTogglesValid
********************/
int CheckTogglesValid(MHProposal *MHp, Network *nwp) {
  int fvalid;
  int i;
  DegreeBound *bd=MHp->bd;

  if(!bd) return 1;

  /* *** don't forget when getting attributes that tail-> head */
  int *tailattr = (int *) Calloc(bd->attrcount, int);
  int *headattr = (int *) Calloc(bd->attrcount, int);
  
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
  
  Free(tailattr);
  Free(headattr);
  
  /* Undo proposed toggles (of edges(tail, head)) */
  for (i=0; i<MHp->ntoggles; i++)
    ToggleEdge(MHp->toggletail[i], MHp->togglehead[i], nwp);

  return fvalid;
}

int CheckConstrainedTogglesValid(MHProposal *MHp, Network *nwp)
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
    int *tailattr = (int *) Calloc(bd->attrcount, int);
    int *headattr = (int *) Calloc(bd->attrcount, int);
    
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
    Free(tailattr);
    Free(headattr);
  }
  /* Make proposed toggles (of edges (tail, head), not (head, tail) */
  for (i=0; i<MHp->ntoggles; i++)
    ToggleEdge(MHp->toggletail[i], MHp->togglehead[i], nwp);
  
  return fvalid;
}
