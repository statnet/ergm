/*  File src/MHproposals_block.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2013 Statnet Commons
 */
#include "MHproposals.h"
#include "ergm_edgelist.h"
#include "ergm_changestat.h"
#include "ergm_MHblockdiag.h"
#include "ergm_MHstorage.h"

/* Shorthand. */

/*********************
 void MH_blockdiag

 Block-diagonal sampling
*********************/
MH_I_FN(Mi_blockdiag){
  ALLOC_STORAGE(1, MH_BlockDiagSampInfo, b);
  double *inputs = MH_INPUTS; // Need an throw-away variable since unpacker updates the position. 
  *b = unpack_BlockDiagSampInfo(&inputs, BIPARTITE, DIRECTED);
  MHp->ntoggles=1;
}

MH_P_FN(Mp_blockdiag){
  GET_STORAGE(MH_BlockDiagSampInfo, b);

  BD_LOOP({
      GetRandDyadBlockDiag(Mtail, Mhead, b);
    });
}

/********************
   void MH_blockTNT

   Block-diagonal TNT sampling
***********************/
MH_I_FN(Mi_blockdiagTNT){
  ALLOC_STORAGE(1, MH_BlockDiagSampInfo, b);
  double *inputs = MH_INPUTS+1; // Need an throw-away variable since unpacker updates the position. 
  *b = unpack_BlockDiagSampInfo(&inputs, BIPARTITE, DIRECTED);
  MHp->ntoggles=1;
}

MH_P_FN(Mp_blockdiagTNT){
  GET_STORAGE(MH_BlockDiagSampInfo, b);

  const double comp=0.5, odds = comp/(1.0-comp);
  
  Dyad ndyads = MH_INPUTS[0];
  Edge nedges=nwp->nedges;
  
  double logratio=0; 

  BD_LOOP({
      if (unif_rand() < comp && nedges > 0) { /* Select a tie at random */
	// Note that, by construction, this tie will be within a block.
	GetRandEdge(Mtail, Mhead, nwp);
	/* Thanks to Robert Goudie for pointing out an error in the previous 
	   version of this sampler when proposing to go from nedges==0 to nedges==1 
	   or vice versa.  Note that this happens extremely rarely unless the 
	   network is small or the parameter values lead to extremely sparse 
	   networks.  */
	logratio = log((nedges==1 ? 1.0/(comp*ndyads + (1.0-comp)) :
			 nedges / (odds*ndyads + nedges)));
      }else{ /* Select a dyad at random within a block */
	GetRandDyadBlockDiag(Mtail, Mhead, b);
	
	if(IS_OUTEDGE(Mtail[0],Mhead[0])!=0){
	  logratio = log((nedges==1 ? 1.0/(comp*ndyads + (1.0-comp)) :
				nedges / (odds*ndyads + nedges)));
	}else{
	  logratio = log((nedges==0 ? comp*ndyads + (1.0-comp) :
				1.0 + (odds*ndyads)/(nedges + 1)));
	}
      }
    });
  MHp->logratio += logratio;
}
