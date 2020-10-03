/*  File src/MHproposals_block.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2020 Statnet Commons
 */
#include "MHproposals.h"
#include "ergm_edgelist.h"

/*********************
 void MH_blockdiag

 Block-diagonal sampling
*********************/
MH_P_FN(MH_blockdiag){  

  /* *** don't forget tail-> head now */

  Vertex tail, head;
  static Vertex blks;
  static double *blkpos, *blkcwt; 
  
  if(MHp->ntoggles == 0) { /* Initialize randomtoggle */
    blks=MHp->inputs[0];
    blkpos = MHp->inputs+1;
    blkcwt = MHp->inputs+1+blks+1; 
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

/*********************
 void MH_blockdiagB

 Block-diagonal sampling for bipartite graphs
*********************/
MH_P_FN(MH_blockdiagB){  

  /* *** don't forget tail-> head now */

  static Vertex blks;
  static double *eblkpos, *ablkpos, *blkcwt; 
  
  if(MHp->ntoggles == 0) { /* Initialize randomtoggle */
    blks=MHp->inputs[0];
    eblkpos = MHp->inputs+1;
    ablkpos = MHp->inputs+1+blks+1;
    blkcwt = MHp->inputs+1+blks+1+blks+1;
    MHp->ntoggles=1;
    return;
  }
  
  BD_LOOP({
      double r = unif_rand();
      // TODO: Use bisection to perform this search in O(log b) instead of O(b) time. 
      Vertex blk = 1;
      while(r>blkcwt[blk-1]) blk++;
      Mtail[0] = eblkpos[blk-1]+1 + unif_rand() * (eblkpos[blk]-eblkpos[blk-1]);
      Mhead[0] = ablkpos[blk-1]+1 + unif_rand() * (ablkpos[blk]-ablkpos[blk-1]);
    });
}


/********************
   void MH_blockTNT

   Block-diagonal TNT sampling
***********************/
MH_P_FN(MH_blockdiagTNT)
{
  /* *** don't forget tail-> head now */
  
  Vertex tail, head, blks=MHp->inputs[1];
  double *blkpos = MHp->inputs+2, *blkcwt = MHp->inputs+2+blks+1, logratio=0; 
  Edge nedges=EDGECOUNT(nwp);
  static double P=0.5;
  static double Q, DP, DO;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    Dyad ndyads = MHp->inputs[0];
    Q = 1-P;
    DP = P*ndyads;
    DO = DP/Q;
    return;
  }
  
  BD_LOOP({
      if (unif_rand() < P && nedges > 0) { /* Select a tie at random */
	// Note that, by construction, this tie will be within a block.
	GetRandEdge(Mtail, Mhead, nwp);
	/* Thanks to Robert Goudie for pointing out an error in the previous 
	   version of this sampler when proposing to go from nedges==0 to nedges==1 
	   or vice versa.  Note that this happens extremely rarely unless the 
	   network is small or the parameter values lead to extremely sparse 
	   networks.  */
	logratio = TNT_LR_E(nedges, Q, DP, DO);
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
	  logratio = TNT_LR_DE(nedges, Q, DP, DO);
	}else{
	  logratio = TNT_LR_DN(nedges, Q, DP, DO);
	}
      }
    });
  MHp->logratio += logratio;
}

/********************
   void MH_blockTNTB

   Block-diagonal TNT sampling for bipartite graphs
***********************/
MH_P_FN(MH_blockdiagTNTB)
{
  /* *** don't forget tail-> head now */

  static double *eblkpos, *ablkpos, *blkcwt; 
  static Vertex blks;
  double logratio=0; 
  Edge nedges=EDGECOUNT(nwp);
  static double P=0.5;
  static double Q, DP, DO;

  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    blks=MHp->inputs[1];
    eblkpos = MHp->inputs+2;
    ablkpos = MHp->inputs+2+blks+1;
    blkcwt = MHp->inputs+2+blks+1+blks+1;
    
    Dyad ndyads = MHp->inputs[0];
    Q = 1-P;
    DP = P*ndyads;
    DO = DP/Q;
    return;
  }
  
  BD_LOOP({
      if (unif_rand() < P && nedges > 0) { /* Select a tie at random */
	// Note that, by construction, this tie will be within a block.
	GetRandEdge(Mtail, Mhead, nwp);
	/* Thanks to Robert Goudie for pointing out an error in the previous 
	   version of this sampler when proposing to go from nedges==0 to nedges==1 
	   or vice versa.  Note that this happens extremely rarely unless the 
	   network is small or the parameter values lead to extremely sparse 
	   networks.  */
	logratio = TNT_LR_E(nedges, Q, DP, DO);
      }else{ /* Select a dyad at random within a block */
	double r = unif_rand();
	// TODO: Use bisection to perform this search in O(log b) instead of O(b) time. 
	Vertex blk = 1;
	while(r>blkcwt[blk-1]) blk++;
	Mtail[0] = eblkpos[blk-1]+1 + unif_rand() * (eblkpos[blk]-eblkpos[blk-1]);
	Mhead[0] = ablkpos[blk-1]+1 + unif_rand() * (ablkpos[blk]-ablkpos[blk-1]);

	if(EdgetreeSearch(Mtail[0],Mhead[0],nwp->outedges)!=0){
	  logratio = TNT_LR_DE(nedges, Q, DP, DO);
	}else{
	  logratio = TNT_LR_DN(nedges, Q, DP, DO);
	}
      }
    });
  MHp->logratio += logratio;
}
