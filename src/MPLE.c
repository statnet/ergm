#include "MPLE.h"

/* *****************
 void MPLE_wrapper

 Wrapper for a call from R.  
 Only find the MPLE, so no MCMC stuff is necessary. 

 Re-Written by David Schruth to add compression:
   Since there are often many different 0->1 dyad changes that
   result in identical change statistic vectors, this routine
   now returns only the unique change statistic vectors in the 
   covmat matrix, along with a separate vector, weightsvector,
   giving the number of repetitions of each unique set of change
   statistics.  Note that two sets of statistics are considered
   unique in this context if the corresponding response values
   (i.e., dyad values, edge or no edge) are unequal.
   The value maxNumDyadTypes is the largest allowable number of
   unique sets of change statistics.
 *****************/

void MPLE_wrapper (int *heads, int *tails, int *dnedges,
		   int *dn, int *dflag, int *bipartite, int *nterms, 
		   char **funnames, char **sonames, double *inputs,  
		   int *responsevec, double *covmat,
		   int *weightsvector,
		   double * offset, double * compressedOffset,
		   int maxNumDyadTypes)
{
  Network nw[2];
  Vertex n_nodes = (Vertex) *dn; 
  Edge n_edges = (Edge) *dnedges;
  int directed_flag = *dflag;
  int hammingterm;
  Vertex bip = (Vertex) *bipartite;
  Vertex hhead, htail;
  Edge  nddyads, kedge;
  Model *m;
  ModelTerm *thisterm;

  GetRNGstate(); /* Necessary for R random number generator */
  nw[0]=NetworkInitialize(heads, tails, n_edges, n_nodes, directed_flag, bip, 0);
  m=ModelInitialize(*funnames, *sonames, inputs, *nterms);
  
  hammingterm=ModelTermHamming (*funnames, *nterms);
  if(hammingterm>0){
   Network nwhamming;
   thisterm = m->termarray + hammingterm - 1;
   nddyads = (Edge)(thisterm->inputparams[0]);
   nwhamming=NetworkInitializeD(thisterm->inputparams+1, 
				thisterm->inputparams+1+nddyads,
			       	nddyads, n_nodes, directed_flag, bip,0);
   nddyads=0;
   nw[1]=NetworkInitializeD(thisterm->inputparams+1, 
			   thisterm->inputparams+1+nddyads, nddyads, n_nodes, directed_flag, bip,0);
/*	     Rprintf("made hw[1]\n"); */
   for (kedge=1; kedge <= nwhamming.nedges; kedge++) {
     FindithEdge(&hhead, &htail, kedge, &nwhamming);
     if(EdgetreeSearch(hhead, htail, nw[0].outedges) == 0){
/*	     Rprintf(" in g0 not g hhead %d htail %d\n",hhead, htail); */
       ToggleEdge(hhead, htail, &nw[1]);
     }
   }
   for (kedge=1; kedge <= nw[0].nedges; kedge++) {
     FindithEdge(&hhead, &htail, kedge, &nw[0]);
     if(EdgetreeSearch(hhead, htail, nwhamming.outedges) == 0){
/*	     Rprintf("not g0  in g hhead %d htail %d\n",hhead, htail); */
       ToggleEdge(hhead, htail, &nw[1]);
     }
   }
/*   Rprintf("Initial number of discordant %d Number of g0 ties %d Number of ties in g %d\n",nw[1].nedges, nwhamming.nedges,nw[0].nedges); */
   NetworkDestroy(&nwhamming);
  }

  MpleInitialize(bip, responsevec, covmat, weightsvector,
		 offset, compressedOffset, maxNumDyadTypes, nw, m); 
  
  ModelDestroy(m);
  NetworkDestroy(nw);
  PutRNGstate(); /* Must be called after GetRNGstate before returning to R */
}

int rowsAreSame(double *rowA, double *rowB, int rowLength)
/* service routine for findCovMatRow() compares rowLength elements 
   of two rows for equality.  Return 0 if the rows are not 
   the same, return 1 otherwise. */
{ 
  int i;
  for (i=0; i<rowLength; i++) {
    if(rowA[i]!=rowB[i]) return(0);
  }
  return(1); 
}

int findCovMatRow(double *newRow, double *matrix, int rowLength, int numRows,
	          int *responsevec,
		  double * offset, double * compressedOffset, int curDyadNum)
{
  int r;
  for (r=0; r < numRows; r++) {
    if(rowsAreSame(newRow,matrix+(rowLength*r),rowLength) 
       && (*(responsevec+r)==responsevec[numRows]) 
       && (*(compressedOffset+r)==offset[curDyadNum])){ 
      return(r); }
  }
  return(-1); /* returns -1 if it couldn't find the row (new row is unique) */
}

/*****************
 void MpleInitialize

 For finding the MPLE, an extra bit of initialization is required:  
 we must build the matrix of covariates to be used in the logistic 
 regression routine.  This matrix has #rows equal to the number of 
 possible edges in the network (n choose 2 for an undirected network) 
 and #cols equal to the number of network statistics in the model, 
 which is also the number of parameters.  The row for each edge in 
 the matrix should contain the changes in each of the statistics that 
 would occur if the given edge is toggled from 0 to 1, leaving all 
 other edges as they are in the observed network.  The response vector 
 for the logistic regression is simply the vector of indicators 
 giving the states of the edges in the observed network.
*****************/
void MpleInitialize (Vertex bipartite, int *responsevec, double *covmat, 
		     int *weightsvector,
		     double * offset, double * compressedOffset,
		     int maxNumDyadTypes, Network *nwp, Model *m) {

  int l, d, outflag = 0, inflag = 0, thisRowNumber, thisOffsetNumber,
    foundRowPosition, totalStats, *currentResponse;
  double *thisPreviousRow, *thisCurrentRow, *covMatPosition;
  int curDyadNum;
  Vertex i, j , rowmax;
  ModelTerm *mtp;
  
  covMatPosition = covmat;
  currentResponse = responsevec;
  thisPreviousRow  = thisCurrentRow =  
    (double*) R_alloc(m->n_stats,sizeof(double));
  curDyadNum=0;
  thisRowNumber = 0;
  thisOffsetNumber = 0;
  if(bipartite > 0){
   rowmax=bipartite+1;
  }else{
   rowmax=nwp->nnodes;
  }
  for(i=1; i < rowmax; i++){
    for(j = MAX(i,bipartite)+1; j <= nwp->nnodes; j++){
      for(d=0; d <= nwp->directed_flag; d++){ /*trivial loop if undirected*/
	if (d==1){*currentResponse = inflag = 
		    (EdgetreeSearch(i, j, nwp->inedges) != 0);}
	else{*currentResponse = outflag = 
	       (EdgetreeSearch(i, j, nwp->outedges) != 0);}
	totalStats = 0;
	/* Let mtp loop through each model term */
	for (mtp=m->termarray; mtp < m->termarray + m->n_terms; mtp++){
	  mtp->dstats = covMatPosition + totalStats;
	  /* Now call d_xxx function, which updates mtp->dstats to reflect
	     changing the current dyad.  */
	  if(d==0){
	    (*(mtp->func))(1, &i, &j, mtp, nwp);
	  }
	  else{ 
	    (*(mtp->func))(1, &j, &i, mtp, nwp);
	  }
	  /* dstats values reflect changes in current dyad; for MPLE, 
	     values must reflect going from 0 to 1.  Thus, we have to reverse 
	     the sign of dstats whenever the current edge exists. */
	  if((d==0 && outflag) || (d==1 && inflag)){
	    for(l=0; l<mtp->nstats; l++){
	      mtp->dstats[l] = -mtp->dstats[l];
	    }
	  }
	  /* Update mtp->dstats pointer to skip ahead by mtp->nstats */
	  totalStats += mtp->nstats; 
	}
	/* Check to see if statistics found at covMatPosition match
           any rows already in covmat matrix (along with corresponding
           response values and offset values) */
	foundRowPosition =  findCovMatRow(covMatPosition, covmat, 
					  m->n_stats, thisRowNumber, 
					  responsevec, offset, 
					  compressedOffset, curDyadNum);
	if(foundRowPosition>=0){  /* Not unique */
	  weightsvector[foundRowPosition]++;
	}else{                    /* unique */
	  if(thisRowNumber<maxNumDyadTypes){ 
	    weightsvector[thisRowNumber]=1;
	    compressedOffset[thisRowNumber] = offset[curDyadNum];
	    /* Shift the pointer n parameters forward in
               the covariate matrix vector */
	    covMatPosition += m->n_stats; /* New row in covmat matrix */
	    currentResponse++; /* New response value */
	    thisRowNumber++; /* New # unique rows */
	  } else{ /* Do nothing for now if thisRowNumber >=maxNumDyadTypes */ }
	}
	curDyadNum++;
      }
    }
  }
}



