/*  File src/MHproposals.c in package ergm, part of the Statnet suite
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
#include "ergm_changestat.h"
#include "ergm_rlebdm.h"
#include "ergm_MHstorage.h"
#include "ergm_unsorted_edgelist.h"
#include "ergm_weighted_population.h"
#include "ergm_dyadgen.h"
#include "ergm_Rutil.h"
#include "ergm_BDStrat_proposals.h"
#include "ergm_hash_edgelist.h"
#include "ergm_nodelist.h"

/*********************
 void MH_randomtoggle

 Default MH algorithm
*********************/
MH_P_FN(MH_randomtoggle){  

  /* *** don't forget tail-> head now */

  if(MHp->ntoggles == 0) { /* Initialize randomtoggle */
    MH_STORAGE = DyadGenInitializeR(MHp->R, nwp, FALSE);
    MHp->ntoggles=1;
    return;
  }
  
  BD_LOOP({
      DyadGenRandDyad(Mtail, Mhead, MH_STORAGE);
    });
}

MH_F_FN(Mf_randomtoggle){
  DyadGenDestroy(MH_STORAGE);
  MH_STORAGE = NULL;
}

/********************
   void MH_TNT
   Propose ONLY edges on a static list
   Use TNT weights.
   This is a fusion of MH_DissolutionMLETNT and MH_TNT:

   A "intersect" network is constructed that is the intersection of
   dyads on the static list and the edges present in nwp. Then,
   standard TNT procedure is followed, but the dyad space (and the
   number of dyads) is the number of dyads in the static list and the
   network for the ties is the ties in the discord network.
***********************/

MH_P_FN(Mp_TNT){
  if(MHp->ntoggles == 0) { /* Initialize randomtoggle */
    MH_STORAGE = DyadGenInitializeR(MHp->R, nwp, TRUE);
    MHp->ntoggles=1;
    return;
  }

  const double P=0.5, Q=1-P;
  double DP = P*((DyadGen *)MH_STORAGE)->ndyads, DO = DP/Q;

  Edge nedges = DyadGenEdgecount(MH_STORAGE);
  double logratio=0;
  BD_LOOP({
      if (unif_rand() < P && nedges > 0) { /* Select a tie at random from the network of eligibles */
        DyadGenRandEdge(Mtail, Mhead, MH_STORAGE);
	logratio = TNT_LR_E(nedges, Q, DP, DO);
      }else{ /* Select a dyad at random from the list */
	DyadGenRandDyad(Mtail, Mhead, MH_STORAGE);
	
	if(IS_OUTEDGE(Mtail[0],Mhead[0])){
	  logratio = TNT_LR_DE(nedges, Q, DP, DO);
	}else{
	  logratio = TNT_LR_DN(nedges, Q, DP, DO);
	}
      }
    });
  MHp->logratio += logratio;
}

MH_F_FN(Mf_TNT){
  DyadGenDestroy(MH_STORAGE);
  MH_STORAGE = NULL;
}


/********************
    MH_BDStratTNT
********************/

// struct definition in ergm_BDStrat_proposals.h

MH_I_FN(Mi_BDStratTNT) {
  // process the inputs and initialize all the edgelists in storage; set MHp->ntoggles to 1
  MHp->ntoggles = 1;

  ALLOC_STORAGE(1, BDStratTNTStorage, sto);

  sto->CD = getListElement(getListElement(MHp->R, "flags"), "CD") != R_NilValue;
  
  sto->bound = asInteger(getListElement(MHp->R, "bound"));
  sto->nmixtypes = asInteger(getListElement(MHp->R, "nmixtypes"));
  sto->nstratlevels = asInteger(getListElement(MHp->R, "nattrcodes"));
  int nbdlevels = asInteger(getListElement(MHp->R, "bd_levels"));

  sto->mixtypestoupdate = Calloc(sto->nmixtypes, int);
  
  sto->strat_vattr = INTEGER(getListElement(MHp->R, "strat_vattr"));

  sto->bd_vattr = INTEGER(getListElement(MHp->R, "bd_vattr"));
    
  sto->nodelist = NodeListInitialize(sto->bound, 
                                     sto->strat_vattr, 
                                     sto->nstratlevels, 
                                     sto->nmixtypes, 
                                     INTEGER(getListElement(MHp->R, "strattailattrs")), 
                                     INTEGER(getListElement(MHp->R, "stratheadattrs")), 
                                     NULL, 
                                     sto->bd_vattr, 
                                     nbdlevels, 
                                     INTEGER(getListElement(MHp->R, "bd_mixtypes"))[0], 
                                     INTEGER(getListElement(MHp->R, "bd_mixtypes"))[1], 
                                     INTEGER(getListElement(MHp->R, "bd_tails")), 
                                     INTEGER(getListElement(MHp->R, "bd_heads")), 
                                     NULL,
                                     INTEGER(getListElement(MHp->R, "nodecountsbypairedcode")), 
                                     nwp);
    
  sto->strat_vattr--; // so node indices line up correctly
  sto->bd_vattr--; // so node indices line up correctly  
    
  UnsrtEL **els = Calloc(sto->nmixtypes, UnsrtEL *);
  for(int i = 0; i < sto->nmixtypes; i++) {
    els[i] = UnsrtELInitialize(0, NULL, NULL, FALSE);
  }
      
  sto->indmat = Calloc(sto->nstratlevels, int *);
  sto->indmat[0] = INTEGER(getListElement(MHp->R, "indmat"));
  for(int i = 1; i < sto->nstratlevels; i++) {
    sto->indmat[i] = sto->indmat[i - 1] + sto->nstratlevels;
  }

  int **amat = Calloc(nbdlevels, int *);
  amat[0] = INTEGER(getListElement(MHp->R, "amat"));
  for(int i = 1; i < nbdlevels; i++) {
    amat[i] = amat[i - 1] + nbdlevels;
  }
  
  EXEC_THROUGH_NET_EDGES(tail, head, e, {
    int index = sto->indmat[sto->strat_vattr[tail]][sto->strat_vattr[head]];
    int allowed = amat[sto->bd_vattr[tail]][sto->bd_vattr[head]];
    if(index >= 0 && allowed) {
      UnsrtELInsert(tail, head, els[index]);
    }  
  });
  Free(amat);

  sto->hash = Calloc(sto->nmixtypes, HashEL *);
  for(int i = 0; i < sto->nmixtypes; i++) {
    sto->hash[i] = HashELInitialize(els[i]->nedges, els[i]->tails + 1, els[i]->heads + 1, FALSE, DIRECTED);
  }

  sto->originalprobvec = Calloc(sto->nmixtypes, double);
  int empirical_flag = asInteger(getListElement(MHp->R, "empirical_flag"));
  if(empirical_flag) {
    for(int i = 0; i < sto->nmixtypes; i++) {
      if(sto->hash[i]->list->nedges > 0) {
        sto->originalprobvec[i] = sto->hash[i]->list->nedges;
      }
    }
  } else {
    memcpy(sto->originalprobvec, REAL(getListElement(MHp->R, "probvec")), sto->nmixtypes*sizeof(double));
  }

  // determine what mixing types are initially toggleable 
  double *currentprobvec = Calloc(sto->nmixtypes, double);  
  for(int i = 0; i < sto->nmixtypes; i++) {
    // if any edges or dyads of this type are toggleable
    if(sto->hash[i]->list->nedges > 0 || NodeListDyadCountPositive(sto->nodelist, i)) {
      currentprobvec[i] = sto->originalprobvec[i];
      sto->currentcumprob += sto->originalprobvec[i];
    }
  }
    
  sto->wtp = WtPopInitialize(sto->nmixtypes, currentprobvec);
  Free(currentprobvec);

  // zero proposal probability is an error
  if(WtPopSumWts(sto->wtp) == 0) {
    MHp->ntoggles = MH_FAILED;
    return;
  }

  for(Vertex vertex = 1; vertex <= N_NODES; vertex++) {
    if(IN_DEG[vertex] + OUT_DEG[vertex] > sto->bound) {
      error("degree bound is violated by initial network; proposal cannot proceed");
    }
  }
}

MH_P_FN(MH_BDStratTNT) {
  GET_STORAGE(BDStratTNTStorage, sto);

  // sample a toggleable strat mixing type on which to make a proposal
  sto->stratmixingtype = WtPopGetRand(sto->wtp);
  
  // number of edges of this mixing type
  int nedgestype = sto->hash[sto->stratmixingtype]->list->nedges;
  
  Dyad ndyadstype = NodeListDyadCount(sto->nodelist, sto->stratmixingtype);
  
  int edgeflag;
  
  if((unif_rand() < 0.5 && nedgestype > 0) || ndyadstype == 0) {
    // propose toggling off an existing edge of strat mixing type sto->stratmixingtype
    HashELGetRand(Mtail, Mhead, sto->hash[sto->stratmixingtype]);    
    edgeflag = TRUE;
  } else {
    // select a random BD toggleable dyad of strat mixing type sto->stratmixingtype and propose toggling it
    NodeListGetRandWithCount(Mtail, Mhead, sto->nodelist, sto->stratmixingtype, ndyadstype);

    // now check if the dyad we drew is already an edge or not
    edgeflag = IS_OUTEDGE(Mtail[0],Mhead[0]);
  }
  
  sto->tailmaxl = IN_DEG[Mtail[0]] + OUT_DEG[Mtail[0]] == sto->bound - 1 + edgeflag;
  sto->headmaxl = IN_DEG[Mhead[0]] + OUT_DEG[Mhead[0]] == sto->bound - 1 + edgeflag;
  
  // compute proposed dyad count for current mixing type (only)
  Dyad proposeddyadstype = NodeListDyadCountOnToggle(*Mtail, *Mhead, sto->nodelist, sto->stratmixingtype, edgeflag ? 1 : -1, sto->tailmaxl, sto->headmaxl);
  
  ComputeChangesToToggleability(Mtail, Mhead, edgeflag, sto);
  
  double prob_weight = sto->currentcumprob/sto->proposedcumprob;
  
  // the rationale for the logratio is similar to that given for BDTNT, with two additional considerations:
  // 
  // - counts of edges and BD toggleable dyads should be for the current mixing type only, and
  //
  // - it is possible for a strat mixing type to reach zero toggleable dyads (by having no edges and also no BD toggleable dyads);
  //   such a strat mixing type cannot be selected when we choose the strat mixing type at the top of the P_FN, and so we must 
  //   "disable" it until it comes to have toggleable dyads again; the term prob_weight adjusts for the fact that the total weight 
  //   given to toggleable strat mixing types may be different in the current and proposed networks, and thus the probability
  //   to select the current strat mixing type may be different in the current and proposed networks
  
  if(edgeflag) {
    MHp->logratio = log(prob_weight*(((nedgestype == 1 ? 1.0 : 0.5)/proposeddyadstype))/(((ndyadstype == 0 ? 1.0/nedgestype : (0.5/nedgestype) + (sto->tailmaxl || sto->headmaxl ? 0.0 : 0.5/ndyadstype)))));
  } else {
    MHp->logratio = log(prob_weight*((proposeddyadstype == 0 ? 1.0/(nedgestype + 1) : (0.5/(nedgestype + 1)) + (sto->tailmaxl || sto->headmaxl ? 0.0 : 0.5/proposeddyadstype))/((nedgestype == 0 ? 1.0 : 0.5)/ndyadstype)));
  }
}

MH_U_FN(Mu_BDStratTNT) {
  GET_STORAGE(BDStratTNTStorage, sto);
  if(sto->CD) {
    sto->stratmixingtype = sto->indmat[sto->strat_vattr[tail]][sto->strat_vattr[head]];

    sto->tailmaxl = IN_DEG[tail] + OUT_DEG[tail] == sto->bound - 1 + edgeflag;
    sto->headmaxl = IN_DEG[head] + OUT_DEG[head] == sto->bound - 1 + edgeflag;

    ComputeChangesToToggleability(&tail, &head, edgeflag, sto);    
  }
  // update edgelist
  HashELToggleKnown(tail, head, sto->hash[sto->stratmixingtype], edgeflag);

  // update nodelists as needed
  NodeListToggleKnownIf(tail, head, sto->nodelist, !edgeflag, sto->tailmaxl, sto->headmaxl);
  
  // if any strat mixing types have changed toggleability status, update prob info accordingly
  if(sto->nmixtypestoupdate > 0) {
    sto->currentcumprob = sto->proposedcumprob;
    for(int i = 0; i < sto->nmixtypestoupdate; i++) {
      WtPopSetWt(sto->mixtypestoupdate[i], edgeflag ? sto->originalprobvec[sto->mixtypestoupdate[i]] : 0, sto->wtp);          
    }
  }
}

MH_F_FN(Mf_BDStratTNT) {
  // Free all the things
  GET_STORAGE(BDStratTNTStorage, sto);

  NodeListDestroy(sto->nodelist);    
  
  for(int i = 0; i < sto->nmixtypes; i++) {
    HashELDestroy(sto->hash[i]);
  }
  Free(sto->hash);

  Free(sto->indmat);

  Free(sto->originalprobvec);

  Free(sto->mixtypestoupdate);  

  WtPopDestroy(sto->wtp);
  // MHp->storage itself should be Freed by MHProposalDestroy
}

/********************
    MH_StratTNT
********************/

// struct definition in ergm_BDStrat_proposals.h

MH_I_FN(Mi_StratTNT) {
  // process the inputs and initialize all the edgelists in storage; set MHp->ntoggles to 1
  MHp->ntoggles = 1;

  // used during initialization only; not retained in storage
  int nlevels = asInteger(getListElement(MHp->R, "nlevels"));

  ALLOC_STORAGE(1, StratTNTStorage, sto);

  sto->vattr = INTEGER(getListElement(MHp->R, "nodecov"));
  
  sto->CD = getListElement(getListElement(MHp->R, "flags"), "CD") != R_NilValue;
  
  sto->nmixtypes = asInteger(getListElement(MHp->R, "nmixtypes"));
  
  sto->nodelist = NodeListInitialize(N_NODES - 1, // i.e. no bound on degree
                                     sto->vattr, 
                                     nlevels, 
                                     sto->nmixtypes, 
                                     INTEGER(getListElement(MHp->R, "tailattrs")), 
                                     INTEGER(getListElement(MHp->R, "headattrs")), 
                                     INTEGER(getListElement(MHp->R, "nodecountsbycode")),
                                     NULL, 0, 0, 0, NULL, NULL, NULL, NULL, nwp);
  
  sto->els = Calloc(sto->nmixtypes, UnsrtEL *);
  for(int i = 0; i < sto->nmixtypes; i++) {
    sto->els[i] = UnsrtELInitialize(0, NULL, NULL, FALSE);
  }
  
  sto->indmat = Calloc(nlevels, int *);
  sto->indmat[0] = INTEGER(getListElement(MHp->R, "indmat"));
  for(int i = 1; i < nlevels; i++) {
    sto->indmat[i] = sto->indmat[i - 1] + nlevels;
  }
  
  sto->vattr--; // decrement so node indices line up correctly
  EXEC_THROUGH_NET_EDGES(tail, head, e, {
    int index = sto->indmat[sto->vattr[tail]][sto->vattr[head]];
    if(index >= 0) {
      UnsrtELInsert(tail, head, sto->els[index]);
    }      
  });
    
  int empirical_flag = asInteger(getListElement(MHp->R, "empirical"));
  if(empirical_flag) {
    double *probvec = Calloc(sto->nmixtypes, double);
    for(int i = 0; i < sto->nmixtypes; i++) {
      probvec[i] = sto->els[i]->nedges;
    }
    sto->wtp = WtPopInitialize(sto->nmixtypes, probvec);
    Free(probvec);
  } else {
    sto->wtp = WtPopInitialize(sto->nmixtypes, REAL(getListElement(MHp->R, "probvec")));
  }

  // zero total proposal probability is an error
  if(WtPopSumWts(sto->wtp) == 0) {
    MHp->ntoggles = MH_FAILED;
    return;
  }
  
  sto->ndyadstype = Calloc(sto->nmixtypes, Dyad);
  for(int i = 0; i < sto->nmixtypes; i++) {
    sto->ndyadstype[i] = NodeListDyadCount(sto->nodelist, i);
    // positive proposal probability with zero dyads is an error;
    // note that we may wish to relax this condition
    if(sto->ndyadstype[i] == 0 && WtPopGetWt(i, sto->wtp) > 0) {
      MHp->ntoggles = MH_FAILED;
      return;        
    }
  }
}

MH_P_FN(MH_StratTNT) {
  GET_STORAGE(StratTNTStorage, sto);
  
  // sample a strat mixing type on which to make a proposal
  sto->currentmixingtype = WtPopGetRand(sto->wtp);
  
  // number of edges of this mixing type
  int nedgestype = sto->els[sto->currentmixingtype]->nedges;
  
  // number of dyads of this mixing type
  Dyad ndyadstype = sto->ndyadstype[sto->currentmixingtype];
  
  BD_LOOP({
    if(unif_rand() < 0.5 && nedgestype > 0) {
      // select an existing edge of type sto->currentmixingtype at random, and propose toggling it off
      UnsrtELGetRand(Mtail, Mhead, sto->els[sto->currentmixingtype]);
      
      // logratio is essentially copied from TNT, because the probability of 
      // choosing this particular mixing type cancels upon taking the ratio;
      // still need to count only edges and dyads of the appropriate mixing type, though
      MHp->logratio = log((nedgestype == 1 ? 1.0/(0.5*ndyadstype + 0.5) :
                           nedgestype / ((double) ndyadstype + nedgestype)));
    } else {
      // select a dyad of type sto->currentmixingtype and propose toggling it
      NodeListGetRandWithCount(Mtail, Mhead, sto->nodelist, sto->currentmixingtype, ndyadstype);

      if(IS_OUTEDGE(Mtail[0],Mhead[0])) {
        // pick a new edge from the edgelist uniformly at random so we know its index
        // and hence don't have to look up the index of the edge tail -> head; this gives
        // the same probability of picking each existing edge as if we used the tail -> head
        // edge, but also allows us to keep the edgelists unsorted (at the cost of generating
        // an extra random index in this case)
        UnsrtELGetRand(Mtail, Mhead, sto->els[sto->currentmixingtype]);

        MHp->logratio = log((nedgestype == 1 ? 1.0/(0.5*ndyadstype + 0.5) :
                             nedgestype / ((double) ndyadstype + nedgestype)));
      }else{                        
        MHp->logratio = log((nedgestype == 0 ? 0.5*ndyadstype + 0.5 :
                             1.0 + (ndyadstype)/((double) nedgestype + 1)));
      }
    }
  });
}

MH_U_FN(Mu_StratTNT) {
  // add or remove edge from appropriate edgelist
  GET_STORAGE(StratTNTStorage, sto);
  if(sto->CD) {
    sto->currentmixingtype = sto->indmat[sto->vattr[tail]][sto->vattr[head]];
  }
  UnsrtELToggleKnown(tail, head, sto->els[sto->currentmixingtype], edgeflag);
}

MH_F_FN(Mf_StratTNT) {
  // Free all the things
  GET_STORAGE(StratTNTStorage, sto);

  for(int i = 0; i < sto->nmixtypes; i++) {
    UnsrtELDestroy(sto->els[i]);
  }

  NodeListDestroy(sto->nodelist);

  Free(sto->els);
  Free(sto->ndyadstype);
  Free(sto->indmat);
  
  WtPopDestroy(sto->wtp);
  // MHp->storage itself should be Freed by MHProposalDestroy
}

/********************
   void MH_TNT10
   Attempts to do 10 TNT steps at once, but this seems flawed currently
   because it does not correctly update network quantities like nedges
   after each of the 10 proposed toggles.
***********************/
MH_P_FN(MH_TNT10)
{
  /* *** don't forget tail-> head now */
  
  Edge nedges=EDGECOUNT(nwp);
  static double P=0.5;
  static double Q, DP, DO;

  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=10;
    Q = 1-P;
    DP = P*DYADCOUNT(nwp);
    DO = DP/Q;
    return;
  }

  double logratio = 0;
  BD_LOOP({
      logratio = 0;
      for(unsigned int n = 0; n < 10; n++){
	if (unif_rand() < P && nedges > 0) { /* Select a tie at random */
	  GetRandEdge(Mtail, Mhead, nwp);
	  logratio += TNT_LR_E(nedges, Q, DP, DO);
	}else{ /* Select a dyad at random */
	  GetRandDyad(Mtail+n, Mhead+n, nwp);
	  if(IS_OUTEDGE(Mtail[n],Mhead[n])!=0){
	    logratio += TNT_LR_DE(nedges, Q, DP, DO);
	  }else{
	    logratio += TNT_LR_DN(nedges, Q, DP, DO);
	  }
	} 
      }
    });
  MHp->logratio += logratio;
}

/*********************
 void MH_constantedges
 propose pairs of toggles that keep number of edges
 the same.  This is done by (a) choosing an existing edge
 at random; (b) repeatedly choosing dyads at random until 
 one is found that does not have an edge; and (c) proposing
 toggling both these dyads.  Note that step (b) will be very 
 inefficient if the network is nearly complete, so this proposal is
 NOT recommended for such networks.  However, most network
 datasets are sparse, so this is not likely to be an issue.
*********************/
MH_P_FN(MH_ConstantEdges){  
  /* *** don't forget tail-> head now */
  
  if(MHp->ntoggles == 0) { /* Initialize */
    if(nwp->nedges==0 || nwp->nedges==DYADCOUNT(nwp)) MHp->ntoggles=MH_FAILED; /* Empty or full network. */
    else MHp->ntoggles=2;
    return;
  }
  /* Note:  This proposal cannot be used for full or empty observed graphs.
     If desired, we could check for this at initialization phase. 
     (For now, however, no way to easily return an error message and stop.)*/
  BD_LOOP({
      /* First, select edge at random */
      GetRandEdge(Mtail, Mhead, nwp);
      /* Second, select non-edge at random */
      GetRandNonedge(Mtail+1, Mhead+1, nwp);
    });
}

/*********************
 void MH_CondDegreeDist
 It used to be called  MH_CondDegDistSwapToggles
*********************/
MH_P_FN(MH_CondDegreeDist){  
  int noutedge=0, ninedge=0, k, fvalid;
  int k0, j0, j1, k1;
  int j0h, j1h;
  int trynode;
  Vertex e, alter, tail=0, head, head1;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=2;    
    return;
  }

  fvalid = 0;
  trynode = 0;
  while(fvalid==0 && trynode < 500){

  trynode++;
  /* select a node at random */
  while(noutedge+ninedge==0){
    /* select a node at random */
    tail = 1 + unif_rand() * N_NODES;
    ninedge  = nwp->indegree[tail];
    noutedge = nwp->outdegree[tail];
  }

  /* choose a edge of the node at random */
    /* *** don't forget tail-> head now */

  k0 = (int)(unif_rand() * (noutedge+ninedge)); 
  if (k0 < noutedge){
    k=0;
    for(e = EdgetreeMinimum(nwp->outedges, tail);
    ((head = nwp->outedges[e].value) != 0 && k<k0);
    e = EdgetreeSuccessor(nwp->outedges, e)){++k;}
  }else{
    k=0;
    for(e = EdgetreeMinimum(nwp->inedges, tail);
    ((head = nwp->inedges[e].value) != 0 && k<(k0-noutedge));
    e = EdgetreeSuccessor(nwp->inedges, e)){++k;}
  }

  if ( (!DIRECTED && tail > head) ||
  (DIRECTED && k0 >= noutedge) ) {
    Mtail[0] = head;
    Mhead[0] = tail;
  }else{
    Mtail[0] = tail;
    Mhead[0] = head;
  }
  
  k1=0;
  fvalid=0;
  while(fvalid==0 && k1 < 100){
    while((alter = 1 + unif_rand() * N_NODES) == tail);
    fvalid=1;
    if(alter == head){fvalid=0;}
    if (k0 < noutedge || !DIRECTED){
      for(e = EdgetreeMinimum(nwp->outedges, tail);
      (fvalid==1 && ((head1 = nwp->outedges[e].value) != 0));
      e = EdgetreeSuccessor(nwp->outedges, e)){
        if(alter==head1){fvalid=0;}}
    }
    if (k0 >= noutedge || !DIRECTED){
      for(e = EdgetreeMinimum(nwp->inedges, tail);
      (fvalid==1 && ((head1 = nwp->inedges[e].value) != 0));
      e = EdgetreeSuccessor(nwp->inedges, e)){
        if(alter==head1){fvalid=0;}}
    }
    k1++;
  }

  if (k1 == 100){
    fvalid=0;
    continue;
  }
  
  if ( (!DIRECTED && alter > tail) ||
       (DIRECTED && k0 < noutedge) )
    {
      Mtail[1] = tail;
      Mhead[1] = alter;
    }else{
      Mtail[1] = alter;
      Mhead[1] = tail;
    }
  
  if (!DIRECTED){
    /* Check undirected degrees */
    k0 =nwp->outdegree[tail]  + nwp->indegree[tail];
    j0h=nwp->outdegree[head]  + nwp->indegree[head];
    j1h=nwp->outdegree[alter] + nwp->indegree[alter];
    
    j0=j0h-1;
    j1=j1h+1;
    
    if( ( (j0==j1h) && (j1==j0h) ) ){
      fvalid = 1;
    }else{
      fvalid = 0;
    }
  }else{
    /* Check directed degrees */
   if(k0 < noutedge){
     /* Check indegrees */
     j0h=nwp->indegree[head];
     j1h=nwp->indegree[alter];
   }else{
     /* Check outdegrees */
     j0h=nwp->outdegree[head];
     j1h=nwp->outdegree[alter];
   }
   j0=j0h-1;
   j1=j1h+1;
   
   if( ( (j0==j1h) && (j1==j0h) ) ){
     fvalid = 1;
   }else{
     fvalid = 0;
   }
  }
  
  }

  if (trynode==500){
    Mtail[1] = Mtail[0];
    Mhead[1] = Mhead[0];
  }
}

/*********************
 void MH_CondOutDegreeDist
*********************/
MH_P_FN(MH_CondOutDegreeDist){  
  int noutedge=0, k, fvalid=0;
  int k0, k1;
  int trynode;
  Vertex e, alter, tail=0, head, head1;
  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=2;    
    return;
  }

  fvalid = 0;
  trynode = 0;
  while(fvalid==0 && trynode < 1500){

  trynode++;

  while(noutedge==0){
    /* select a node at random */
    tail = 1 + unif_rand() * N_NODES;
    noutedge = nwp->outdegree[tail];
  }
  
  k0 = (int)(unif_rand() * noutedge); 
  k=0;
  for(e = EdgetreeMinimum(nwp->outedges, tail);
      ((head = nwp->outedges[e].value) != 0 && k<k0);
      e = EdgetreeSuccessor(nwp->outedges, e)){++k;}
  Mtail[0] = tail;
  Mhead[0] = head;
  
  k1=0;
  fvalid=0;
  while(fvalid==0 && k1 < 100){
    while((alter = 1 + unif_rand() * N_NODES) == tail);
    fvalid=1;
    if(alter == head){fvalid=0;}
    for(e = EdgetreeMinimum(nwp->outedges, tail);
	(fvalid==1 && ((head1 = nwp->outedges[e].value) != 0));
	e = EdgetreeSuccessor(nwp->outedges, e)){
      if(alter==head1){fvalid=0;}}
    k1++;
  }
  if (k1 == 100){
    fvalid=0;
    continue;
  }
  
  Mtail[1] = tail;
  Mhead[1] = alter;
  }
  
  if(trynode==1500 || !CheckTogglesValid(MHp, nwp)){
      Mtail[0] = 1;
      Mhead[0] = 2;
      Mtail[1] = 1;
      Mhead[1] = 2;
  }
  

}

/*********************
 void MH_CondInDegreeDist
*********************/
MH_P_FN(MH_CondInDegreeDist){  
  int ninedge=0, k, fvalid=0;
  int k0, k1;
  int trynode;
  Vertex e, alter, tail=0, head, head1;

  /* *** don't forget tail-> head now */

  
  if(MHp->ntoggles == 0) { /* Initialize */
    MHp->ntoggles=2;    
    return;
  }

  fvalid = 0;
  trynode = 0;
  while(fvalid==0 && trynode < 1500){

  trynode++;

  while(ninedge==0){
    /* select a node at random */
    tail = 1 + unif_rand() * N_NODES;
    ninedge = nwp->indegree[tail];
  }
  
  k0 = (int)(unif_rand() * ninedge); 
  k=0;
  for(e = EdgetreeMinimum(nwp->inedges, tail);
      ((head = nwp->inedges[e].value) != 0 && k<k0);
      e = EdgetreeSuccessor(nwp->inedges, e)){++k;}
  Mtail[0] = head;
  Mhead[0] = tail;
  
  k1=0;
  fvalid=0;
  while(fvalid==0 && k1 < 100){
    while((alter = 1 + unif_rand() * N_NODES) == tail);
    fvalid=1;
    if(alter == head){fvalid=0;}
    for(e = EdgetreeMinimum(nwp->inedges, tail);
	(fvalid==1 && ((head1 = nwp->inedges[e].value) != 0));
	e = EdgetreeSuccessor(nwp->inedges, e)){
      if(alter==head1){fvalid=0;}}
    k1++;
  }
  if (k1 == 100){
    fvalid=0;
    continue;
  }
  
  Mtail[1] = alter;
  Mhead[1] = tail;
  
  }
  
  if(trynode==1500){
      Mtail[0] = 1;
      Mhead[0] = 2;
      Mtail[1] = 1;
      Mhead[1] = 2;
  }
}

/*********************
 void MH_TwoRandomToggles
*********************/
MH_P_FN(MH_TwoRandomToggles){  
  Vertex tail, head;
  int i;

  /* *** don't forget tail-> head now */
  
  if(MHp->ntoggles == 0) { /* Initialize OneRandomToggle */
    MHp->ntoggles=2;
    return;
  }

  for (i = 0; i < 2; i++){
   tail = 1 + unif_rand() * N_NODES;
   while ((head = 1 + unif_rand() * N_NODES) == tail);
   if (!DIRECTED && tail > head) {
     Mtail[i] = head;
     Mhead[i] = tail;
   }else{
     Mtail[i] = tail;
     Mhead[i] = head;
   }
  }
}

/*********************
 void MH_RandomNode
*********************/
MH_P_FN(MH_randomnode){
  
  Vertex root, alter;
  int j;
  
  if(MHp->ntoggles == 0) { /* Initialize OneRandomToggle */
    MHp->ntoggles= N_NODES - 1;
    return;
  }

  root = 1 + unif_rand() * N_NODES;
  
  j = 0;
  for (alter = 1; alter <= N_NODES; alter++)
    {
      /* there is never an edge (root, root) */
      if (alter != root) {
       if (!DIRECTED && root > alter) {
        Mtail[j] = alter;
        Mhead[j] = root;
       }else{
        Mtail[j] = root;
        Mhead[j] = alter;
       }
       j++;
      }
    }
}

/* The ones below have not been tested */

/*********************
 void MH_ConstrainedCondOutDegDist
*********************/
MH_P_FN(MH_ConstrainedCondOutDegDist){  
  int noutedge=0, k, fvalid=0;
  int k0, k1;
  Vertex e, alter, tail, head, head1;

  /* *** don't forget tail-> head now */

  while(noutedge==0){
    /* select a node at random */
    tail = 1 + unif_rand() * N_NODES;
    noutedge = nwp->outdegree[tail];
  }
  
  k0 = (int)(unif_rand() * noutedge); 
  k=0;
  for(e = EdgetreeMinimum(nwp->outedges, tail);
      ((head = nwp->outedges[e].value) != 0 && k<k0);
      e = EdgetreeSuccessor(nwp->outedges, e)){++k;}
  Mtail[0] = tail;
  Mhead[0] = head;
  
  k1=0;
  fvalid=0;
  while(fvalid==0 && k1 < 100){
    while((alter = 1 + unif_rand() * N_NODES) == tail);
    fvalid=1;
    if(alter == head){fvalid=0;}
    for(e = EdgetreeMinimum(nwp->outedges, tail);
	(fvalid==1 && ((head1 = nwp->outedges[e].value) != 0));
	e = EdgetreeSuccessor(nwp->outedges, e)){
      if(alter==head1){fvalid=0;}}
    k1++;
  }
  if (k1 == 100){
    Mtail[0] = Mhead[0] = 0;
    Mtail[1] = Mhead[1] = 0;
  }
  
  Mtail[1] = tail;
  Mhead[1] = alter;
  
  if (!fvalid){
    Mtail[0] = Mhead[0] = 0;
    Mtail[1] = Mhead[1] = 0;
  }
  
  for(k=0; k < 2; k++){
    if (dEdgeListSearch(Mtail[k], Mhead[k], MH_INPUTS)==0){
      Mtail[0] = Mhead[0] = 0;
      Mtail[1] = Mhead[1] = 0;
    }
  }
}


MH_P_FN(MH_NodePairedTiesToggles){  
  /* chooses a node and toggles all ties and
	 and toggles an equal number of matching nonties
	 for that node */
  int nedge=0,j,k;
  int fvalid = 1;
  Vertex e, tail, prop;

  /* *** don't forget tail-> head now */
  
  /* double to integer coercion */
  tail = 1 + unif_rand() * N_NODES; 
  
  for(e = EdgetreeMinimum(nwp->outedges, tail);
      (prop = nwp->outedges[e].value) != 0; /* loop if */
      e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
    {
      Mtail[nedge] = tail;
      Mhead[nedge] = prop;
      ++nedge;
    }
  for(e = EdgetreeMinimum(nwp->inedges, tail);
      (prop = nwp->inedges[e].value) != 0; /* loop if */
      e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail */
    {
      Mhead[nedge] = tail;
      Mtail[nedge] = prop;
      ++nedge;
    }
  
  if(nedge > N_NODES-nedge){
    Mtail[0] = Mhead[0] = 0;
    Mtail[1] = Mhead[1] = 0;
  }  
  j = 0;
  while (j <=nedge)
    {
      prop = 1 + unif_rand() * N_NODES; 
      k=0;
      fvalid=1;
      while(fvalid==1 && k<nedge+j){
	if(IS_OUTEDGE( MIN(prop,Mtail[k]),
		       MAX(prop,Mtail[k])) +
	   IS_OUTEDGE( MIN(prop,Mhead[k]),
		       MAX(prop,Mhead[k]))==0
	   ){++k;
	}else{
	  fvalid=0;
	}
      }
      if(prop>tail){
	Mtail[j+nedge] = tail;
	Mhead[j+nedge] = prop;
      }else{
	Mtail[j+nedge] = prop;
	Mhead[j+nedge] = tail;
      }
      ++j;
    }
  
  j = 2*nedge;
  if (!CheckTogglesValid(MHp, nwp))
    {
      *Mtail = *Mhead = 0;
    }
}

/*********************
 void MH_OneRandomTnTNode
*********************/
MH_P_FN(MH_OneRandomTnTNode){  
  Vertex tail=0, head, e, head1;
  int noutedge=0, ninedge=0, k0=0, fvalid=0, k;
  /* int ndyad; */

  /* *** don't forget tail-> head now */
  
  /* if ( DIRECTED )
    {
      ndyad = (N_NODES - 1) * N_NODES;
    }else{
      ndyad = (N_NODES - 1) * N_NODES / 2;
    } */

  double logratio=0;
  fvalid=0;
  while(fvalid==0){
    
    if ( unif_rand() < 0.5 && EDGECOUNT(nwp) > 0) 
      {
	
	/* select a tie */
	ninedge=0;
	noutedge=0;
	while(noutedge+ninedge==0){
	  /* select a node at random */
	  tail = 1 + unif_rand() * N_NODES;
	  ninedge = nwp->indegree[tail];
	  noutedge = nwp->outdegree[tail];
	}
	
	k0 = (int)(unif_rand() * (noutedge+ninedge)); 
	if (k0 < noutedge){
	  k=0;
	  for(e = EdgetreeMinimum(nwp->outedges, tail);
	      ((head = nwp->outedges[e].value) != 0 && k<k0);
	      e = EdgetreeSuccessor(nwp->outedges, e)){++k;}
	}else{
	  k=0;
	  for(e = EdgetreeMinimum(nwp->inedges, tail);
	      ((head = nwp->inedges[e].value) != 0 && k<(k0-noutedge));
	      e = EdgetreeSuccessor(nwp->inedges, e)){++k;}
	}
	if ( (!DIRECTED && tail > head) ||
	     (DIRECTED && k0 >= noutedge) )
	  {
	    Mtail[0] = head;
	    Mhead[0] = tail;
	  }else{
	    Mtail[0] = tail;
	    Mhead[0] = head;
	  }
	
	logratio = log(((noutedge+ninedge)*1.0)/(N_NODES-1-noutedge-ninedge-1));
	fvalid =1;
      }else{
	/* Choose random non-tie */

	/* select a node at random */
	ninedge=N_NODES-1;
	noutedge=0;
	while(noutedge+ninedge>=(N_NODES-1)){
	  ninedge=0;
	  /* select a node at random */
	  tail = 1 + unif_rand() * N_NODES;
	  ninedge = nwp->indegree[tail];
	  noutedge = nwp->outdegree[tail];
	}
	
	fvalid=0;
	while(fvalid==0){
	  while ((head = 1 + unif_rand() * N_NODES) == tail);
	  fvalid=1;
	  for(e = EdgetreeMinimum(nwp->outedges, tail);
	      (fvalid==1 && ((head1 = nwp->outedges[e].value) != 0));
	      e = EdgetreeSuccessor(nwp->outedges, e)){
	    if(head==head1){fvalid=0;}}
	  if (!(DIRECTED)){
	    for(e = EdgetreeMinimum(nwp->inedges, tail);
		(fvalid==1 && ((head1 = nwp->inedges[e].value) != 0));
		e = EdgetreeSuccessor(nwp->inedges, e)){
	      if(head==head1){fvalid=0;}}
	  }
	}
	
	if ( (!DIRECTED && tail > head) ||
	     (DIRECTED && k0 >= noutedge) )
	  {
	    Mtail[0] = head;
	    Mhead[0] = tail;
	  }else{
	    Mtail[0] = tail;
	    Mhead[0] = head;
	  }
	
        if ( DIRECTED )
	  {
	    logratio = log((N_NODES-1-noutedge-ninedge)/(noutedge+ninedge+1.0));
	  }else{
	    logratio = log((N_NODES-1-noutedge-ninedge)/(noutedge+ninedge+1.0));
	  }
      }
  }
  MHp->logratio += logratio;
}

/*********************
 void MH_ReallocateWithReplacement
*********************/
MH_P_FN(MH_ReallocateWithReplacement){  
  int i;
  Vertex root;
  Vertex* edges;
  int edgecount = 0;
  
  /* select a node at random */
  root = 1 + unif_rand() * N_NODES;

  edges = (Vertex *) Calloc(N_NODES+1, Vertex);
  for (i = 0; i <= N_NODES; i++)
    edges[i] = NO_EDGE;
  
  /* count current edges and mark them in an array */
  for (i = 1; i <= N_NODES; i++)
    {
      if (root == i) continue;
      if (IS_OUTEDGE(root, i) > 0)
	{
	  edges[i] = OLD_EDGE;
	  edgecount++;
	}
      if (!DIRECTED && (root > i) &&
	  (IS_OUTEDGE(i, root) > 0))
	{
	  edges[i] = OLD_EDGE;
	  edgecount++;
	}
    }
  
  /* select edgecount edges to create */
  for (i = 0; i < edgecount; i++)
    {
      Vertex newhead;
      /* get a new edge, neither the root nor something already chosen */
      while ((newhead = 1 + unif_rand() * N_NODES) == root ||
	     (edges[newhead] & NEW_EDGE))
	;
      
      /* if this edge already exists - (OLD_EDGE | NEW_EDGE) == CAN_IGNORE */
      edges[newhead] = edges[newhead] | NEW_EDGE;
    }
  
  /* index into Mtail/Mhead is  */
  edgecount = 0;
  
  /* add to toggle list:  anything that is non zero in edges array
     should be toggled, whether on or off. */
  for (i = 0; i <= N_NODES; i++)
    {
      if (edges[i] == NO_EDGE || edges[i] == CAN_IGNORE) continue;
      
      /* double to integer coercion */
      Mtail[edgecount] = root;
      Mhead[edgecount] = i;
      
      if (!DIRECTED && (Mtail[edgecount] > Mhead[edgecount]))
	{
	  Vertex temp;
	  temp = Mtail[edgecount];
	  Mtail[edgecount] = Mhead[edgecount];
	  Mhead[edgecount] = temp;
	}
      edgecount++;
    }
  Free(edges);
}

/*********************
 void MH_AllTogglesForOneNode
*********************/
MH_P_FN(MH_AllTogglesForOneNode){
  
  int i;
  int j;
  int root;
  
  root = 1 + unif_rand() * N_NODES;
  
  j = 0;
  for (i = 1; i <= N_NODES; i++)
    {
      /* probability here only do this with .8? */
      
      /* there is never an edge (root, root) */
      if (i == root)
	continue;
      
      /* double to integer coercion */
      Mtail[j] = root;
      Mhead[j] = i;
      
      if (!DIRECTED && (Mtail[j] > Mhead[j]))
	{
	  Vertex temp;
	  temp = Mtail[j];
	  Mtail[j] = Mhead[j];
	  Mhead[j] = temp;
	}
      j++;
    }
}


/*********************
 void MH_SwitchLabelTwoNodesToggles
*********************/
MH_P_FN(MH_SwitchLabelTwoNodesToggles){  
  int nedge1=0, nedge2=0, k, ntoggles;
  Vertex *edges1, *edges2;
  Vertex e, tail2, head2, tail1, head1;

  /* *** don't forget tail-> head now */
  
  /* select a node at random */
  edges1 = (Vertex *) Calloc(N_NODES+1, Vertex);
  edges2 = (Vertex *) Calloc(N_NODES+1, Vertex);
  
  while(nedge1==0){
    tail1 = 1 + unif_rand() * N_NODES;
    
    for(e = EdgetreeMinimum(nwp->outedges, tail1);
	(head1 = nwp->outedges[e].value) != 0; /* loop if */
	e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail1 */
      {
        edges1[nedge1] = head1;
	++nedge1;
      }
    for(e = EdgetreeMinimum(nwp->inedges, tail1);
	(head1 = nwp->inedges[e].value) != 0; /* loop if */
	e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail1 */
      {
        edges1[nedge1] = head1;
	++nedge1;
      }
  }
  
  while((tail2 = 1 + unif_rand() * N_NODES) == tail1);
  
  for(e = EdgetreeMinimum(nwp->outedges, tail2);
      (head2 = nwp->outedges[e].value) != 0; /* loop if */
      e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail2 */
    {
      edges2[nedge2] = head2;
      ++nedge2;
    }
  for(e = EdgetreeMinimum(nwp->inedges, tail2);
      (head2 = nwp->inedges[e].value) != 0; /* loop if */
      e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail2 */
    {
      edges2[nedge2] = head2;
      ++nedge2;
    }
  
  ntoggles = 0;
  for(k=0; k < nedge1; k++){
    if (tail1 > edges1[k])
      {
	Mtail[ntoggles] = edges1[k];
	Mhead[ntoggles] = tail1;
      }
    if (tail1 < edges1[k]){
      Mtail[ntoggles] = tail1;
      Mhead[ntoggles] = edges1[k];
    }
    if(tail1 != edges1[k]) ntoggles++;
  }
  
  for(k=0; k < nedge2; k++){
    if (tail1 > edges2[k])
      {
	Mtail[ntoggles] = edges2[k];
	Mhead[ntoggles] = tail1;
      }
    if (tail1 < edges2[k]){
      Mtail[ntoggles] = tail1;
      Mhead[ntoggles] = edges2[k];
    }
    if(tail1 != edges2[k]) ntoggles++;
  }
  
  for(k=0; k < nedge2; k++){
    if (tail2 > edges2[k])
      {
	Mtail[ntoggles] = edges2[k];
	Mhead[ntoggles] = tail2;
      }
    if (tail2 < edges2[k]){
      Mtail[ntoggles] = tail2;
      Mhead[ntoggles] = edges2[k];
    }
    if(tail2 != edges2[k]) ntoggles++;
  }
  
  for(k=0; k < nedge1; k++){
    if (tail2 > edges1[k])
      {
	Mtail[ntoggles] = edges1[k];
	Mhead[ntoggles] = tail2;
      }
    if (tail2 < edges1[k]){
      Mtail[ntoggles] = tail2;
      Mhead[ntoggles] = edges1[k];
    }
    if(tail2 != edges1[k]) ntoggles++;
  }
  Free(edges1);
  Free(edges2);
}


/*********************
 void MH_ConstrainedCondDegDist
*********************/
MH_P_FN(MH_ConstrainedCondDegDist){  
  int noutedge=0, ninedge=0, k, fvalid=0;
  int k0, j0, j1, k1;
  int j0h, j1h;
  Vertex *outedges, *inedges;
  Vertex e, alter, tail=0, head;

  /* *** don't forget tail-> head now */
  
  /* select a node at random */
  outedges = (Vertex *) Calloc(N_NODES+1, Vertex);
  inedges = (Vertex *) Calloc(N_NODES+1, Vertex);
  
  while(noutedge==0 && ninedge==0){
    tail = 1 + unif_rand() * N_NODES;
    
    for(e = EdgetreeMinimum(nwp->outedges, tail);
	(head = nwp->outedges[e].value) != 0; /* loop if */
	e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
      {
        outedges[noutedge] = head;
	++noutedge;
      }
    for(e = EdgetreeMinimum(nwp->inedges, tail);
	(head = nwp->inedges[e].value) != 0; /* loop if */
	e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail */
      {
        inedges[ninedge] = head;
	++ninedge;
      }
  }
  
  k0 = (int)(unif_rand() * (noutedge+ninedge)); 
  if (k0 < noutedge){
    head = outedges[k0]; 
  }else{
    head = inedges[k0-noutedge]; 
  }
  if ( (!DIRECTED && tail > head) ||
       (  DIRECTED  && k0 >= noutedge) )
    {
      Mtail[0] = head;
      Mhead[0] = tail;
    }else{
      Mtail[0] = tail;
      Mhead[0] = head;
    }
  
  if (dEdgeListSearch(Mtail[0], Mhead[0], MH_INPUTS)==0){
    Mtail[0] = Mhead[0] = 0;
    Mtail[1] = Mhead[1] = 0;
  }
  
  fvalid=0;
  k1=0;
  while(fvalid==0 && k1 < 100){
    while((alter = 1 + unif_rand() * N_NODES) == tail);
    if(alter != head){fvalid=1;}
    fvalid=1;
    if (k0 < noutedge || !(DIRECTED)){
      k=0;
      while(fvalid==1 && noutedge > 0 && k <= noutedge-1){
	if(alter == outedges[k]){fvalid=0;}else{++k;}
      }
    }
    if (k0 >= noutedge || !(DIRECTED)){
      k=0;
      while(fvalid==1 && ninedge > 0 && k <= ninedge-1){
	if(alter == inedges[k]){fvalid=0;}else{++k;}
      }
    }
    k1++;
  }
  
  if (k1 == 100){
    Mtail[0] = Mhead[0] = 0;
    Mtail[1] = Mhead[1] = 0;
  }
  
  if ( (!DIRECTED && alter > tail) ||
       (DIRECTED && k0 < noutedge) )
    {
      Mtail[1] = tail;
      Mhead[1] = alter;
    }else{
      Mtail[1] = alter;
      Mhead[1] = tail;
    }
  
  if (dEdgeListSearch(Mtail[1], Mhead[1], MH_INPUTS)==0){
    Mtail[0] = Mhead[0] = 0;
    Mtail[1] = Mhead[1] = 0;
  }
  
  Free(outedges);
  Free(inedges);
  
  /* Check undirected degrees */

  /* *** don't forget tail-> head now */

  if (!DIRECTED){
    k0=nwp->outdegree[tail]+ nwp->indegree[tail];
    j0h=nwp->outdegree[head]+ nwp->indegree[head];
    j1h=nwp->outdegree[alter]+ nwp->indegree[alter];
    
    j0=j0h-1;
    j1=j1h+1;
    
    if( ( (j0==j1h) && (j1==j0h) ) ){
      fvalid = 1;
    }else{
      fvalid = 0;
    }
  }else{
    if(k0 < noutedge){
      /* Check indegrees */
      j0h=nwp->indegree[head];
      j1h=nwp->indegree[alter];
    }else{
      /* Check outdegrees */
      j0h=nwp->outdegree[head];
      j1h=nwp->outdegree[alter];
    }
    j0=j0h-1;
    j1=j1h+1;
    
    if( ( (j0==j1h) && (j1==j0h) ) ){
      fvalid = 1;
    }else{
      fvalid = 0;
    }
  }
  
  if (!fvalid){
    Mtail[0] = Mhead[0] = 0;
    Mtail[1] = Mhead[1] = 0;
  }
}

void MH_ConstrainedNodePairedTiesToggles (MHProposal *MHp,
       	 Network *nwp) {  
  /* chooses a node and toggles all ties and
     and toggles an equal number of matching nonties
     for that node */
  int nedge=0,j,k;
  int fvalid = 1;
  Vertex e, tail, prop;

  /* *** don't forget tail-> head now */
  
  /* double to integer coercion */
  tail = 1 + unif_rand() * N_NODES; 
  
  for(e = EdgetreeMinimum(nwp->outedges, tail);
      (prop = nwp->outedges[e].value) != 0; /* loop if */
      e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
    {
      Mtail[nedge] = tail;
      Mhead[nedge] = prop;
      ++nedge;
    }
  for(e = EdgetreeMinimum(nwp->inedges, tail);
      (prop = nwp->inedges[e].value) != 0; /* loop if */
      e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail */
    {
      Mhead[nedge] = tail;
      Mtail[nedge] = prop;
      ++nedge;
    }
  
  if(nedge > N_NODES-nedge){
    Mtail[0] = Mhead[0] = 0;
    Mtail[1] = Mhead[1] = 0;
  }  
  j = 0;
  while (j <=nedge)
    {
      prop = 1 + unif_rand() * N_NODES; 
      k=0;
      fvalid=1;
      while(fvalid==1 && k<nedge+j){
	if(IS_OUTEDGE(MIN(prop,Mtail[k]),
			   MAX(prop,Mtail[k])) +
	   IS_OUTEDGE( MIN(prop,Mhead[k]),
			   MAX(prop,Mhead[k]))==0
	   ){++k;
	}else{
	  fvalid=0;}
      }
      if(prop>tail){
	Mtail[j+nedge] = tail;
	Mhead[j+nedge] = prop;
      }else{
	Mtail[j+nedge] = prop;
	Mhead[j+nedge] = tail;
      }
      ++j;
    }
  
  j = 2*nedge;
  if (!CheckConstrainedTogglesValid(MHp, nwp))
    {
      *Mtail = *Mhead = 0;
    }
}

/*********************
 void MH_ConstrainedReallocateWithReplacement
*********************/
void MH_ConstrainedReallocateWithReplacement (MHProposal *MHp,
       	 Network *nwp) {  
  int i;
  Vertex root;
  Vertex* edges;
  int edgecount = 0;
  
  /* select a node at random */
  root = 1 + unif_rand() * N_NODES;

  edges = (Vertex *) Calloc(N_NODES+1, Vertex);
  for (i = 0; i <= N_NODES; i++)
    edges[i] = NO_EDGE;
  
  /* count current edges and mark them in an array */
  for (i = 1; i <= N_NODES; i++)
    {
      if (root == i) continue;
      if (IS_OUTEDGE(root, i) > 0)
	{
	  edges[i] = OLD_EDGE;
	  edgecount++;
	}
      if (!DIRECTED && (root > i) &&
	  (IS_OUTEDGE(i, root) > 0))
	{
	  edges[i] = OLD_EDGE;
	  edgecount++;
	}
    }
  
  /* select edgecount edges to create */
  for (i = 0; i < edgecount; i++)
    {
      Vertex newhead;
      
      /* get a new edge, neither the root nor something already chosen */
      while ((newhead = 1 + unif_rand() * N_NODES) == root ||
	     (edges[newhead] & NEW_EDGE))
	;
      
      /* if this edge already exists - (OLD_EDGE | NEW_EDGE) == CAN_IGNORE */
      edges[newhead] = edges[newhead] | NEW_EDGE;
    }
  
  /* index into Mtail/Mhead is  */
  edgecount = 0;
  
  /* add to toggle list:  anything that is non zero in edges array
     should be toggled, whether on or off. */
  for (i = 0; i <= N_NODES; i++)
    {
      if (edges[i] == NO_EDGE || edges[i] == CAN_IGNORE) continue;
      
      /* double to integer coercion */
      Mtail[edgecount] = root;
      Mhead[edgecount] = i;
      
      if (!DIRECTED && (Mtail[edgecount] > Mhead[edgecount]))
	{
	  Vertex temp;
	  temp = Mtail[edgecount];
	  Mtail[edgecount] = Mhead[edgecount];
	  Mhead[edgecount] = temp;
	}
      edgecount++;
    }
  Free(edges);
}

/*********************
 void MH_ConstrainedAllTogglesForOneNode
*********************/
void MH_ConstrainedAllTogglesForOneNode (MHProposal *MHp,
					 Network *nwp) {
  int i;
  int j;
  int root;
  
  root = 1 + unif_rand() * N_NODES;
  
  j = 0;
  for (i = 1; i <= N_NODES; i++)
    {
      /* probability here only do this with .8? */
      
      /* there is never an edge (root, root) */
      if (i == root)
	continue;
      
      /* double to integer coercion */
      Mtail[j] = root;
      Mhead[j] = i;
      
      if (!DIRECTED && (Mtail[j] > Mhead[j]))
	{
	  Vertex temp;
	  temp = Mtail[j];
	  Mtail[j] = Mhead[j];
	  Mhead[j] = temp;
	}
      j++;
    }
}

/*********************
 void MH_ConstrainedTwoRandomToggles
*********************/
void MH_ConstrainedTwoRandomToggles (MHProposal *MHp,
				 Network *nwp) {  
  int i;
  
  for (i = 0; i < 2; i++)
    {
      /* double to integer coercion */
      Mtail[i] = 1 + unif_rand() * N_NODES; 
      while ((Mhead[i] = 1 + unif_rand() * N_NODES) == Mtail[i]);
      
      while(dEdgeListSearch(Mtail[i], Mhead[i], MH_INPUTS)==0){
	Mtail[i] = 1 + unif_rand() * N_NODES; 
	while ((Mhead[i] = 1 + unif_rand() * N_NODES) == Mtail[i]);
      }
      if (!DIRECTED && Mtail[i] > Mhead[i]) 
	{
	  Vertex temp;
	  temp = Mtail[i];
	  Mtail[i] = Mhead[i];
	  Mhead[i] = temp;
	}
    }
  
  if (!CheckConstrainedTogglesValid(MHp, nwp))
    {
      Mtail[0] = Mhead[0] = 0;
      Mtail[1] = Mhead[1] = 0;
    }  
}

/*********************
 void MH_ConstrainedCondDeg
*********************/
void MH_ConstrainedCondDeg (MHProposal *MHp,
					 Network *nwp) {  
  /* WARNING: THIS NEEDS TO BE FIXED */
  int nedge1=0, nedge2=0, k, toomany, fvalid=0;
  Vertex *edges1, *edges2;
  Vertex e, tail2=0, head2, tail1, head1;
  
  /* select a node at random */
  edges1 = (Vertex *) Calloc(N_NODES+1, Vertex);
  edges2 = (Vertex *) Calloc(N_NODES+1, Vertex);
  
  while(nedge1==0){
    tail1 = 1 + unif_rand() * N_NODES;
    
    for(e = EdgetreeMinimum(nwp->outedges, tail1);
	(head1 = nwp->outedges[e].value) != 0; /* loop if */
	e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail1 */
      {
        edges1[nedge1] = head1;
	++nedge1;
      }
    for(e = EdgetreeMinimum(nwp->inedges, tail1);
	(head1 = nwp->inedges[e].value) != 0; /* loop if */
	e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail1 */
      {
        edges1[nedge1] = head1;
	++nedge1;
      }
  }
  
  head1 = edges1[(int)(unif_rand() * nedge1)]; 
  if (tail1 > head1)
    {
      Mtail[0] = head1;
      Mhead[0] = tail1;
    }else{
      Mtail[0] = tail1;
      Mhead[0] = head1;
    }
   
  toomany = 0;
  while(nedge2==0 && toomany < 100){
    fvalid=0;
    while(fvalid==0){
      while((tail2 = 1 + unif_rand() * N_NODES) == tail1);
      k=0;
      fvalid=1;
      while(fvalid==1 && k < nedge1){
	if(tail2 == edges1[k]){fvalid=0;}else{++k;}
      }
    }

    for(e = EdgetreeMinimum(nwp->outedges, tail2);
	(head2 = nwp->outedges[e].value) != 0; /* loop if */
	e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail2 */
      {
        edges2[nedge2] = head2;
	++nedge2;
      }
    for(e = EdgetreeMinimum(nwp->inedges, tail2);
	(head2 = nwp->inedges[e].value) != 0; /* loop if */
	e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail2 */
      {
        edges2[nedge2] = head2;
	++nedge2;
      }
    ++toomany;
  }
  if (toomany==100){
    Mtail[0] = Mhead[0] = 0;
    Mtail[1] = Mhead[1] = 0;
  }
  toomany=0;
  fvalid=0;
  while(fvalid==0 && toomany < 10){
    while((head2 = edges2[(int)(unif_rand() * nedge2)]) == tail1);
    k=0;
    fvalid=1;
    while(fvalid==1 && k < nedge1){
      if(head2 == edges1[k]){fvalid=0;}else{++k;}
    }
    ++toomany;
  }
  if (!fvalid || toomany==10){
    Mtail[0] = Mhead[0] = 0;
    Mtail[1] = Mhead[1] = 0;
    Free(edges1);
    Free(edges2);
      }
  if (tail2 > head2)
    {
      Mtail[1] = head2;
      Mhead[1] = tail2;
    }else{
      Mtail[1] = tail2;
      Mhead[1] = head2;
    }
  Free(edges1);
  Free(edges2);
}

/*********************
 void MH_ConstrainedSwitchLabelTwoNodesToggles
*********************/
void MH_ConstrainedSwitchLabelTwoNodesToggles (MHProposal *MHp,
       	 Network *nwp)  {  
  int nedge1=0, nedge2=0, k, ntoggles;
  Vertex *edges1, *edges2;
  Vertex e, tail2, head2, tail1, head1;

  /* *** don't forget tail-> head now */
  
  /* select a node at random */

  edges1 = (Vertex *) Calloc(N_NODES+1, Vertex);
  edges2 = (Vertex *) Calloc(N_NODES+1, Vertex);

  while(nedge1==0){
    tail1 = 1 + unif_rand() * N_NODES;
    
    for(e = EdgetreeMinimum(nwp->outedges, tail1);
	(head1 = nwp->outedges[e].value) != 0; /* loop if */
	e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail1 */
      {
        edges1[nedge1] = head1;
	++nedge1;
      }
    for(e = EdgetreeMinimum(nwp->inedges, tail1);
	(head1 = nwp->inedges[e].value) != 0; /* loop if */
	e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail1 */
      {
        edges1[nedge1] = head1;
	++nedge1;
      }
  }
  
  while((tail2 = 1 + unif_rand() * N_NODES) == tail1);
  
  for(e = EdgetreeMinimum(nwp->outedges, tail2);
      (head2 = nwp->outedges[e].value) != 0; /* loop if */
      e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail2 */
    {
      edges2[nedge2] = head2;
      ++nedge2;
    }
  for(e = EdgetreeMinimum(nwp->inedges, tail2);
      (head2 = nwp->inedges[e].value) != 0; /* loop if */
      e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail2 */
    {
      edges2[nedge2] = head2;
      ++nedge2;
    }
  
  ntoggles = 0;
  for(k=0; k < nedge1; k++){
    if (tail1 > edges1[k])
      {
	Mtail[ntoggles] = edges1[k];
	Mhead[ntoggles] = tail1;
      }
    if (tail1 < edges1[k]){
      Mtail[ntoggles] = tail1;
      Mhead[ntoggles] = edges1[k];
    }
    if(tail1 != edges1[k]) ntoggles++;
  }
  
  for(k=0; k < nedge2; k++){
    if (tail1 > edges2[k])
      {
	Mtail[ntoggles] = edges2[k];
	Mhead[ntoggles] = tail1;
      }
    if (tail1 < edges2[k]){
      Mtail[ntoggles] = tail1;
      Mhead[ntoggles] = edges2[k];
    }
    if(tail1 != edges2[k]) ntoggles++;
  }
  
  for(k=0; k < nedge2; k++){
    if (tail2 > edges2[k])
      {
	Mtail[ntoggles] = edges2[k];
	Mhead[ntoggles] = tail2;
      }
    if (tail2 < edges2[k]){
      Mtail[ntoggles] = tail2;
      Mhead[ntoggles] = edges2[k];
    }
    if(tail2 != edges2[k]) ntoggles++;
  }
  
  for(k=0; k < nedge1; k++){
    if (tail2 > edges1[k])
      {
	Mtail[ntoggles] = edges1[k];
	Mhead[ntoggles] = tail2;
      }
    if (tail2 < edges1[k]){
      Mtail[ntoggles] = tail2;
      Mhead[ntoggles] = edges1[k];
    }
    if(tail2 != edges1[k]) ntoggles++;
  }
  Free(edges1);
  Free(edges2);
}

/*********************
 void MH_ConstantEdgesToggles
*********************/
MH_P_FN(MH_ConstantEdgesToggles){  
  int noutedge=0, ninedge=0, k, fvalid=0;
  int k0, k1;
  Vertex e, alter, tail, head, head1;

  /* *** don't forget tail-> head now */
  
  while(noutedge+ninedge==0){
    /* select a node at random */
    tail = 1 + unif_rand() * N_NODES;
    ninedge  = nwp->indegree[tail];
    noutedge = nwp->outdegree[tail];
  }
  
  k0 = (int)(unif_rand() * (noutedge+ninedge)); 
  if (k0 < noutedge){
    k=0;
    for(e = EdgetreeMinimum(nwp->outedges, tail);
	((head = nwp->outedges[e].value) != 0 && k<k0);
	e = EdgetreeSuccessor(nwp->outedges, e)){++k;}
  }else{
    k=0;
    for(e = EdgetreeMinimum(nwp->inedges, tail);
	((head = nwp->inedges[e].value) != 0 && k<(k0-noutedge));
	e = EdgetreeSuccessor(nwp->inedges, e)){++k;}
  }
  
  if ( (!DIRECTED && tail > head) ||
       (DIRECTED && k0 >= noutedge) )
    {
      Mtail[0] = head;
      Mhead[0] = tail;
    }else{
      Mtail[0] = tail;
      Mhead[0] = head;
    }
  
  k1=0;
  fvalid=0;
  while(fvalid==0 && k1 < 100){
    while((alter = 1 + unif_rand() * N_NODES) == tail);
    fvalid=1;
    if(alter == head){fvalid=0;}
    if (k0 < noutedge || !(DIRECTED)){
      for(e = EdgetreeMinimum(nwp->outedges, tail);
	  (fvalid==1 && ((head1 = nwp->outedges[e].value) != 0));
	  e = EdgetreeSuccessor(nwp->outedges, e)){
	if(alter==head1){fvalid=0;}}
    }
    if (k0 >= noutedge || !(DIRECTED)){
      for(e = EdgetreeMinimum(nwp->inedges, tail);
	  (fvalid==1 && ((head1 = nwp->inedges[e].value) != 0));
	  e = EdgetreeSuccessor(nwp->inedges, e)){
	if(alter==head1){fvalid=0;}}
    }
    k1++;
  }
  if (k1 == 100){
    Mtail[0] = Mhead[0] = 0;
    Mtail[1] = Mhead[1] = 0;
  }
  
  if ( (!DIRECTED && alter > tail) ||
       (DIRECTED && k0 < noutedge) )
    {
      Mtail[1] = tail;
      Mhead[1] = alter;
    }else{
      Mtail[1] = alter;
      Mhead[1] = tail;
    }
  
  if (!fvalid){
    Mtail[0] = Mhead[0] = 0;
    Mtail[1] = Mhead[1] = 0;
  }else{  
  }
}

/*********************
 void MH_CondDegSwitchToggles
*********************/
MH_P_FN(MH_CondDegSwitchToggles){  
  int noutedge, ninedge, i;
  int k, k0, toomany;
  Vertex e, tail, head;

  /* *** don't forget tail-> head now */
  
  /* select a node at random */
  for (i = 0; i < 2; i++){
    toomany=0;
    noutedge=0;
    ninedge=0;
    while(noutedge==0 && ninedge==0 && toomany < 100){
      tail = 1 + unif_rand() * N_NODES;
      ninedge=0;
      noutedge=0;
      while(noutedge+ninedge==0){
	/* select a node at random */
	tail = 1 + unif_rand() * N_NODES;
	ninedge = nwp->indegree[tail];
	noutedge = nwp->outdegree[tail];
      }
      ++toomany;
    }
    
    if (toomany == 100){
      Mtail[0] = Mhead[0] = 0;
      Mtail[1] = Mhead[1] = 0;
    }
    
    k0 = (int)(unif_rand() * (noutedge+ninedge)); 
    if (k0 < noutedge){
      k=0;
      for(e = EdgetreeMinimum(nwp->outedges, tail);
	  ((head = nwp->outedges[e].value) != 0 && k<k0);
	  e = EdgetreeSuccessor(nwp->outedges, e)){++k;}
    }else{
      k=0;
      for(e = EdgetreeMinimum(nwp->inedges, tail);
	  ((head = nwp->inedges[e].value) != 0 && k<(k0-noutedge));
	  e = EdgetreeSuccessor(nwp->inedges, e)){++k;}
    }
    if ( (!DIRECTED && tail > head) ||
	 (DIRECTED && k0 >= noutedge) )
      {
	Mtail[i] = head;
	Mhead[i] = tail;
      }else{
	Mtail[i] = tail;
	Mhead[i] = head;
      }
  }
  
  if (IS_OUTEDGE( Mtail[0],Mhead[1]) ||
      IS_OUTEDGE( Mtail[1],Mhead[0]) ){
    Mtail[0] = Mhead[0] = 0;
    Mtail[1] = Mhead[1] = 0;
  }
  
  if ( (!DIRECTED && Mtail[0] > Mhead[1]) )
    {
      Mtail[2] = Mhead[1];
      Mhead[2] = Mtail[0];
    }else{
      Mtail[2] = Mtail[0];
      Mhead[2] = Mhead[1];
    }
  
  if ( (!DIRECTED && Mtail[1] > Mhead[0]) )
    {
      Mtail[3] = Mhead[0];
      Mhead[3] = Mtail[1];
    }else{
      Mtail[3] = Mtail[1];
      Mhead[3] = Mhead[0];
    }
}


