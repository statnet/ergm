/*  File src/MHproposals_dyadnoise.c in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2023 Statnet Commons
 */
#include "ergm_MHproposal.h"
#include "ergm_MHproposal_bd.h"
#include "ergm_changestat.h"
#include "ergm_edgelist.h"
#include "ergm_MHstorage.h"

/********************
   void MH_dyadnoiseTNT
   Tie/no tie:  Gives at least 50% chance of
   proposing a toggle of an existing edge, as opposed
   to simple random toggles that rarely do so in sparse 
   networks

   Takes an "observed" network and a list of precomputed ratios, and
   tweaks the probability of the network by
   p(nwp[i,j])to(obs[i,j]). (Network is passed second.)

   TODO: Proposals could almost certainly be more efficient.

***********************/
MH_P_FN(MH_dyadnoiseTNT)
{
  /* *** don't forget tail-> head now */
  static double comp=0.5;
  static double odds, o0s0, o0s1, o1s0, o1s1;
  static Edge ndyads;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MH_STORAGE = DegreeBoundInitializeR(MHp->R, nwp);
    MHp->ntoggles=1;
    odds = comp/(1.0-comp);
    ndyads = DYADCOUNT(N_NODES, BIPARTITE, DIRECTED);
    // o=observed, s=state
    o0s0 = MHp->inputs[0]; o0s1 = MHp->inputs[1];
    o1s0 = MHp->inputs[2]; o1s1 = MHp->inputs[3];
    return;
  }
  
  BD_LOOP(MH_STORAGE, {
      if (unif_rand() < comp && N_EDGES > 0) { /* Select a tie at random */
	GetRandEdge(Mtail, Mhead, nwp);
	/* Thanks to Robert Goudie for pointing out an error in the previous 
	   version of this sampler when proposing to go from N_EDGES==0 to N_EDGES==1 
	   or vice versa.  Note that this happens extremely rarely unless the 
	   network is small or the parameter values lead to extremely sparse 
	   networks.  */
	MHp->logratio += log((N_EDGES==1 ? 1.0/(comp*ndyads + (1.0-comp)) :
			      N_EDGES / (odds*ndyads + N_EDGES)));

	int obs = dEdgeListSearch(Mtail[0],Mhead[0],MHp->inputs+4)!=0;	
	MHp->logratio += obs?o1s1:o0s1;
      }else{ /* Select a dyad at random */
	GetRandDyad(Mtail, Mhead, nwp);

	int obs = dEdgeListSearch(Mtail[0],Mhead[0],MHp->inputs+4)!=0;	

	if(IS_OUTEDGE(Mtail[0],Mhead[0])!=0){
	  MHp->logratio += log((N_EDGES==1 ? 1.0/(comp*ndyads + (1.0-comp)) :
				N_EDGES / (odds*ndyads + N_EDGES)));
	  MHp->logratio += obs?o1s1:o0s1;
	}else{
	  MHp->logratio += log((N_EDGES==0 ? comp*ndyads + (1.0-comp) :
				1.0 + (odds*ndyads)/(N_EDGES + 1)));
	  MHp->logratio += obs?o1s0:o0s0;
	}
      }
    });
}

MH_F_FN(Mf_dyadnoiseTNT){
  DegreeBoundDestroy(MH_STORAGE);
}


/********************
   void MH_dyadnoisemTNT
   Tie/no tie:  Gives at least 50% chance of
   proposing a toggle of an existing edge, as opposed
   to simple random toggles that rarely do so in sparse 
   networks

   Takes an "observed" network and a list of precomputed ratio *matrices*, and
   tweaks the probability of the network by
   p(nwp[i,j])to(obs[i,j]). (Network is passed second.)

   TODO: Proposals could almost certainly be more efficient.

***********************/
MH_I_FN(Mi_dyadnoisemTNT){
  MH_STORAGE = DegreeBoundInitializeR(MHp->R, nwp);
  MHp->ntoggles = 1;
}

MH_P_FN(MH_dyadnoisemTNT)
{
  /* *** don't forget tail-> head now */  
  static double comp=0.5;
  static double odds, *o0s0, *o0s1, *o1s0, *o1s1, *onwp;
  static Edge ndyads;
  static Vertex ntails;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    odds = comp/(1.0-comp);
    ndyads = DYADCOUNT(N_NODES, BIPARTITE, DIRECTED);
    // o=observed, s=state
    unsigned int matsize = BIPARTITE? (N_NODES-BIPARTITE)*BIPARTITE : N_NODES*N_NODES;
    o0s0 = MHp->inputs; o0s1 = MHp->inputs+matsize;
    o1s0 = MHp->inputs+matsize*2; o1s1 = MHp->inputs+matsize*3;
    onwp = MHp->inputs+matsize*4;
    ntails = BIPARTITE? BIPARTITE:N_NODES;
    return;
  }
  
  BD_LOOP(MH_STORAGE, {
      if (unif_rand() < comp && N_EDGES > 0) { /* Select a tie at random */
	GetRandEdge(Mtail, Mhead, nwp);
	/* Thanks to Robert Goudie for pointing out an error in the previous 
	   version of this sampler when proposing to go from N_EDGES==0 to N_EDGES==1 
	   or vice versa.  Note that this happens extremely rarely unless the 
	   network is small or the parameter values lead to extremely sparse 
	   networks.  */
	MHp->logratio += log((N_EDGES==1 ? 1.0/(comp*ndyads + (1.0-comp)) :
			      N_EDGES / (odds*ndyads + N_EDGES)));

	int obs = dEdgeListSearch(Mtail[0],Mhead[0],onwp)!=0;	
	MHp->logratio += (obs?o1s1:o0s1)[(Mtail[0]-1) + (Mhead[0]-1)*ntails]; // R matrix serialization is column-major.
      }else{ /* Select a dyad at random */
	GetRandDyad(Mtail, Mhead, nwp);

	int obs = dEdgeListSearch(Mtail[0],Mhead[0],onwp)!=0;	

	if(IS_OUTEDGE(Mtail[0],Mhead[0])!=0){
	  MHp->logratio += log((N_EDGES==1 ? 1.0/(comp*ndyads + (1.0-comp)) :
				N_EDGES / (odds*ndyads + N_EDGES)));
	  MHp->logratio += (obs?o1s1:o0s1)[(Mtail[0]-1) + (Mhead[0]-1)*ntails]; // R matrix serialization is column-major.
	}else{
	  MHp->logratio += log((N_EDGES==0 ? comp*ndyads + (1.0-comp) :
				1.0 + (odds*ndyads)/(N_EDGES + 1)));
	  MHp->logratio += (obs?o1s0:o0s0)[(Mtail[0]-1) + (Mhead[0]-1)*ntails]; // R matrix serialization is column-major.
	}
      }
    });
}

MH_F_FN(Mf_dyadnoisemTNT){
  DegreeBoundDestroy(MH_STORAGE);
}


/********************
   void MH_dyadnoise

   Takes an "observed" network and a list of precomputed ratios, and
   tweaks the probability of the network by
   p(nwp[i,j])to(obs[i,j]). (Network is passed second.)

   TODO: Proposals could almost certainly be more efficient.

***********************/
MH_I_FN(Mi_dyadnoise){
  MH_STORAGE = DegreeBoundInitializeR(MHp->R, nwp);
  MHp->ntoggles = 1;
}

MH_P_FN(MH_dyadnoise)
{
  /* *** don't forget tail-> head now */
  
  static double o0s0, o0s1, o1s0, o1s1;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    // o=observed, s=state
    o0s0 = MHp->inputs[0]; o0s1 = MHp->inputs[1];
    o1s0 = MHp->inputs[2]; o1s1 = MHp->inputs[3];
    return;
  }
  
  BD_LOOP(MH_STORAGE, {
      GetRandDyad(Mtail, Mhead, nwp);
      
      int obs = dEdgeListSearch(Mtail[0],Mhead[0],MHp->inputs+4)!=0;	
      
      if(IS_OUTEDGE(Mtail[0],Mhead[0])!=0){
	MHp->logratio += obs?o1s1:o0s1;
      }else{
	MHp->logratio += obs?o1s0:o0s0;
      }
    });
}

MH_F_FN(Mf_dyadnoise){
  DegreeBoundDestroy(MH_STORAGE);
}


/********************
   void MH_dyadnoisem

   Takes an "observed" network and a list of precomputed ratio *matrices*, and
   tweaks the probability of the network by
   p(nwp[i,j])to(obs[i,j]). (Network is passed second.)

   TODO: Proposals could almost certainly be more efficient.

***********************/
MH_I_FN(Mi_dyadnoisem){
  MH_STORAGE = DegreeBoundInitializeR(MHp->R, nwp);
  MHp->ntoggles = 1;
}

MH_P_FN(MH_dyadnoisem)
{
  /* *** don't forget tail-> head now */
  
  static double *o0s0, *o0s1, *o1s0, *o1s1, *onwp;
  static Vertex ntails;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=1;
    // o=observed, s=state
    unsigned int matsize = BIPARTITE? (N_NODES-BIPARTITE)*BIPARTITE : N_NODES*N_NODES;
    o0s0 = MHp->inputs; o0s1 = MHp->inputs+matsize;
    o1s0 = MHp->inputs+matsize*2; o1s1 = MHp->inputs+matsize*3;
    onwp = MHp->inputs+matsize*4;
    ntails = BIPARTITE? BIPARTITE:N_NODES;
    return;
  }
  
  BD_LOOP(MH_STORAGE, {
      GetRandDyad(Mtail, Mhead, nwp);

      int obs = dEdgeListSearch(Mtail[0],Mhead[0],onwp)!=0;	

      if(IS_OUTEDGE(Mtail[0],Mhead[0])!=0){
	MHp->logratio += (obs?o1s1:o0s1)[(Mtail[0]-1) + (Mhead[0]-1)*ntails]; // R matrix serialization is column-major.
      }else{
	MHp->logratio += (obs?o1s0:o0s0)[(Mtail[0]-1) + (Mhead[0]-1)*ntails]; // R matrix serialization is column-major.
      }
    });
}

MH_F_FN(Mf_dyadnoisem){
  DegreeBoundDestroy(MH_STORAGE);
}
