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
#include "edgelist.h"

/* Shorthand. */
#define Mtail (MHp->toggletail)
#define Mhead (MHp->togglehead)

/*********************
 void MH_blockdiag

 Block-diagonal sampling
*********************/
void MH_blockdiag (MHproposal *MHp, Network *nwp)  {  

  /* *** don't forget tail-> head now */

  Vertex tail, head, blks=MHp->inputs[0];
  double *blkpos = MHp->inputs+1, *blkcwt = MHp->inputs+1+blks+1; 
  
  if(MHp->ntoggles == 0) { /* Initialize randomtoggle */
    MHp->ntoggles=1;
    return;
  }
  
  BD_LOOP({
      double r = unif_rand();
      // TODO: Use bisection to perform this search in O(log b) instead of O(b) time. 
      Vertex blk = 1;
      while(r>blkcwt[blk-1]) blk++;
      tail = blkpos[blk-1]+1 + unif_rand() * (blkpos[blk]-blkpos[blk-1]);
      while ((head = blkpos[blk-1]+1 + unif_rand() * (blkpos[blk]-blkpos[blk-1])) == tail);
      
      if (!nwp->directed_flag && tail > head) {
	Mtail[0] = head;
	Mhead[0] = tail;
      }else{
	Mtail[0] = tail;
	Mhead[0] = head;
      }
    });
}

/********************
   void MH_blockTNT

   Block-diagonal TNT sampling
***********************/
void MH_blockdiagTNT (MHproposal *MHp, Network *nwp) 
{
  /* *** don't forget tail-> head now */
  
  Vertex tail, head, blks=MHp->inputs[1];
  double *blkpos = MHp->inputs+2, *blkcwt = MHp->inputs+2+blks+1, logratio=0; 
  Edge nedges=nwp->nedges;
  static double comp=0.5;
  static double odds;
  static Edge ndyads;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    odds = comp/(1.0-comp);
    ndyads = MHp->inputs[0];
    return;
  }
  
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
	double r = unif_rand();
	// TODO: Use bisection to perform this search in O(log b) instead of O(b) time. 
	Vertex blk = 1;
	while(r>blkcwt[blk-1]) blk++;
	tail = blkpos[blk-1]+1 + unif_rand() * (blkpos[blk]-blkpos[blk-1]);
	while ((head = blkpos[blk-1]+1 + unif_rand() * (blkpos[blk]-blkpos[blk-1])) == tail);
	
	if (tail > head && !nwp->directed_flag)  {
	  Mtail[0] = head;
	  Mhead[0] = tail;
	}else{
	  Mtail[0] = tail;
	  Mhead[0] = head;
	}
	if(EdgetreeSearch(Mtail[0],Mhead[0],nwp->outedges)!=0){
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
