/*  File src/MHproposal_bd.c in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2022 Statnet Commons
 */
#include "ergm_MHproposal_bd.h"
#include "ergm_changestat.h"

/***********************
 DegreeBound* DegreeBoundInitializeR
************************/

DegreeBound* DegreeBoundInitializeR(SEXP MHpR, Network *nwp){
  SEXP bdR = getListElement(MHpR, "bd");
  if(!isNULL(bdR)){
    return(
           DegreeBoundInitialize(INTEGER(getListElement(bdR,"attribs")),
                                 INTEGER(getListElement(bdR,"maxout")),
                                 INTEGER(getListElement(bdR,"maxin")),
                                 INTEGER(getListElement(bdR,"minout")),
                                 INTEGER(getListElement(bdR,"minin")),
                                 asInteger(getListElement(bdR,"condAllDegExact")),
                                 length(getListElement(bdR,"attribs")),
                                 nwp)
           );
  }else return(NULL);
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
  if(!bd) return;
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
int CheckTogglesValid(DegreeBound *bd, MHProposal *MHp, Network *nwp) {
  int fvalid;
  int i;

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

        EXEC_THROUGH_FOUTEDGES(MHp->toggletail[i], e, v, {
	    for (k=0; k < bd->attrcount; k++)
	      if (bd->attribs[v-1 + k*nwp->nnodes]) tailattr[k]++;
	  });

        /* calculate head indegree totals for each attribute
        for each inedge of the head */

	EXEC_THROUGH_FINEDGES(MHp->togglehead[i], e, v, {
	    for (k=0; k < bd->attrcount; k++)
	      if (bd->attribs[v-1 + k*nwp->nnodes]) headattr[k]++;
	  });

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

	EXEC_THROUGH_EDGES(MHp->toggletail[i], e, v, {
	    for (k=0; k < bd->attrcount; k++)
	      if (bd->attribs[v-1 + k*nwp->nnodes])
		tailattr[k]++;
	  });

	      /* calculate head totals for each attribute
        for each outedge and inedge of the head */

	EXEC_THROUGH_EDGES(MHp->togglehead[i], e, v, {
	    for (k=0; k < bd->attrcount; k++)
	      if (bd->attribs[v-1 + k*nwp->nnodes])
		headattr[k]++;
	  });

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

int CheckConstrainedTogglesValid(DegreeBound *bd, MHProposal *MHp, Network *nwp)
{
  int fvalid = 1;
  int i;

  if(!bd) return 1;

  /* Make proposed toggles */
  for (i=0; i<MHp->ntoggles; i++)
    ToggleEdge(MHp->toggletail[i], MHp->togglehead[i], nwp);

  /* if we're bounding degrees by attribute */
  if (bd->fBoundDegByAttr && fvalid)
  {
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

	EXEC_THROUGH_FOUTEDGES(MHp->toggletail[i], e, v, {
	    for (k=0; k < bd->attrcount; k++)
	      if (bd->attribs[v-1 + k*nwp->nnodes])
		tailattr[k]++;
	  });

	      /* calculate head indegree totals for each attribute
        for each inedge of the head */

	EXEC_THROUGH_FINEDGES(MHp->togglehead[i], e, v, {
	    for (k=0; k < bd->attrcount; k++)
	      if (bd->attribs[v-1 + k*nwp->nnodes])
		headattr[k]++;
	  });

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

	EXEC_THROUGH_EDGES(MHp->toggletail[i], e, v, {
	    for (k=0; k < bd->attrcount; k++)
	      if (bd->attribs[v-1 + k*nwp->nnodes])
		tailattr[k]++;
	  });

	      /* calculate head totals for each attribute
        for each outedge and inedge of the head */
	EXEC_THROUGH_EDGES(MHp->togglehead[i], e, v, {
	    for (k=0; k < bd->attrcount; k++)
	      if (bd->attribs[v-1 + k*nwp->nnodes])
		headattr[k]++;
	  });

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
