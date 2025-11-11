/*  File src/MHproposals.c in package ergm, part of the Statnet suite of
 *  packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#include "MHproposals.h"
#include "ergm_edgelist.h"
#include "ergm_changestat.h"
#include "ergm_rlebdm.h"
#include "ergm_MHstorage.h"
#include "ergm_unsorted_edgelist.h"
#include "ergm_weighted_population.h"
#include "ergm_Rutil.h"
#include "ergm_BDStrat_proposals.h"
#include "ergm_hash_edgelist.h"
#include "ergm_BDNodeLists.h"
#include "ergm_BDStratBlocks.h"

/*********************
 void MH_randomtoggle

 Default MH algorithm with dyad generator API.
*********************/
MH_I_FN(Mi_randomtoggle){
  INIT_DYADGEN_AND_DEGREE_BOUND_AND_MODEL(FALSE);
  MHp->ntoggles = storage->gen->ndyads!=0 ? 1 : MH_FAILED;
}

MH_P_FN(MH_randomtoggle){
  GET_STORAGE(StoreDyadGenAndDegreeBoundAndModel, storage);
  DyadGenRandDyad(Mtail, Mhead, storage->gen);
  CHECK_BD(storage->bd);
  CHECK_CHANGESTATS(storage->m);
}

MH_F_FN(Mf_randomtoggle){
  DESTROY_DYADGEN_AND_DEGREE_BOUND_AND_MODEL;
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

MH_I_FN(Mi_TNT){
  INIT_DYADGEN_AND_DEGREE_BOUND_AND_MODEL(TRUE);
  MHp->ntoggles = storage->gen->ndyads!=0 ? 1 : MH_FAILED;
}

MH_P_FN(Mp_TNT){
  GET_STORAGE(StoreDyadGenAndDegreeBoundAndModel, storage);

  const double P=0.5, Q=1-P;
  double DP = P*storage->gen->ndyads, DO = DP/Q;

  Edge nedges = DyadGenEdgecount(storage->gen);
  double logratio=0;
  bool edgestate;
  if (unif_rand() < P && nedges > 0) { /* Select a tie at random from the network of eligibles */
    DyadGenRandEdge(Mtail, Mhead, storage->gen);
    logratio = TNT_LR_E(nedges, Q, DP, DO);
    edgestate = true;
  }else{ /* Select a dyad at random from the list */
    DyadGenRandDyad(Mtail, Mhead, storage->gen);

    if((edgestate = IS_OUTEDGE(Mtail[0],Mhead[0]))){
      logratio = TNT_LR_DE(nedges, Q, DP, DO);
    }else{
      logratio = TNT_LR_DN(nedges, Q, DP, DO);
    }
  }

  CHECK_BD(storage->bd);
  CHECK_CHANGESTATS(storage->m, edgestate);

  MHp->logratio += logratio;
}


MH_F_FN(Mf_TNT){
  DESTROY_DYADGEN_AND_DEGREE_BOUND_AND_MODEL;
}

/********************
    MH_BDStratTNT
********************/

MH_I_FN(Mi_BDStratTNT) {
  // single-toggle proposal
  MHp->ntoggles = 1;

  // BDStratTNTStorage struct definition in inst/ergm_BDStrat_proposals.h
  ALLOC_STORAGE(1, BDStratTNTStorage, sto);

  // store attribute pointers; decrement so nodal indices line up correctly
  sto->strat_vattr = INTEGER(getListElement(MHp->R, "strat_vattr")) - 1;
  sto->blocks_vattr = INTEGER(getListElement(MHp->R, "blocks_vattr")) - 1;
  sto->bd_vattr = INTEGER(getListElement(MHp->R, "bd_vattr")) - 1;

  // read/store number of levels for each attribute
  sto->strat_nlevels = asInteger(getListElement(MHp->R, "strat_nlevels"));
  int nblockslevels = asInteger(getListElement(MHp->R, "blocks_nlevels"));
  sto->bd_nlevels = asInteger(getListElement(MHp->R, "bd_nlevels"));

  // set up degree bound arrays
  sto->maxout = R_Calloc(sto->bd_nlevels, int *);
  sto->maxin = DIRECTED ? R_Calloc(sto->bd_nlevels, int *) : sto->maxout;
  sto->maxout[0] = INTEGER(getListElement(MHp->R, "maxout")) - 1;
  if(DIRECTED) {
    sto->maxin[0] = INTEGER(getListElement(MHp->R, "maxin")) - 1;
  }
  for(int i = 1; i < sto->bd_nlevels; i++) {
    sto->maxout[i] = sto->maxout[i - 1] + N_NODES;
    if(DIRECTED) {
      sto->maxin[i] = sto->maxin[i - 1] + N_NODES;
    }
  }

  // initialize degree storage
  sto->indegree = R_Calloc(sto->bd_nlevels, int *);
  sto->outdegree = R_Calloc(sto->bd_nlevels, int *);
  for(int i = 0; i < sto->bd_nlevels; i++) {
    sto->indegree[i] = R_Calloc(N_NODES + 1, int);
    sto->outdegree[i] = R_Calloc(N_NODES + 1, int);
  }

  // tabulate initial degrees
  EXEC_THROUGH_NET_EDGES(tail, head, e, {
    sto->indegree[sto->bd_vattr[tail]][head]++;
    sto->outdegree[sto->bd_vattr[head]][tail]++;
  });

  // odds and ends
  sto->CD = getListElement(getListElement(MHp->R, "flags"), "CD") != R_NilValue;
  sto->strat_nmixtypes = length(getListElement(MHp->R, "probvec"));
  sto->strat_mixtypestoupdate = R_Calloc(sto->strat_nmixtypes, int);

  // set up submaximal node lists given degrees, degree bounds, and attribute information
  sto->lists = BDNodeListsInitialize(sto->maxout,
                                     sto->maxin,
                                     sto->outdegree,
                                     sto->indegree,
                                     INTEGER(getListElement(MHp->R, "combined_vattr")) - 1,
                                     asInteger(getListElement(MHp->R, "combined_nlevels")),
                                     sto->bd_vattr,
                                     sto->bd_nlevels,
                                     INTEGER(getListElement(MHp->R, "combined_vattr_counts")),
                                     nwp);

  // set up blocks sampler (submaximal node lists form margins of blocks)
  sto->blocks = BDStratBlocksInitialize(sto->lists,
                                        sto->strat_nlevels,
                                        sto->strat_nmixtypes,
                                        INTEGER(getListElement(MHp->R, "strat_tails")),
                                        INTEGER(getListElement(MHp->R, "strat_heads")),
                                        nblockslevels,
                                        INTEGER(getListElement(MHp->R, "blocks_nmixtypes")),
                                        INTEGER(getListElement(MHp->R, "blocks_tails")),
                                        INTEGER(getListElement(MHp->R, "blocks_heads")),
                                        sto->bd_nlevels,
                                        INTEGER(getListElement(MHp->R, "bd_nmixtypes")),
                                        INTEGER(getListElement(MHp->R, "bd_tails")),
                                        INTEGER(getListElement(MHp->R, "bd_heads")),
                                        nwp);

  // set up stratified edge storage
  int *strattailattrs = INTEGER(getListElement(MHp->R, "strat_tails"));
  int *stratheadattrs = INTEGER(getListElement(MHp->R, "strat_heads"));

  // indmat stores indices of included strat mixing types
  sto->indmat = R_Calloc(sto->strat_nlevels, int *);
  for(int i = 0; i < sto->strat_nlevels; i++) {
    sto->indmat[i] = R_Calloc(sto->strat_nlevels, int);
    for(int j = 0; j < sto->strat_nlevels; j++) {
      sto->indmat[i][j] = -1;
    }
  }
  for(int i = 0; i < sto->strat_nmixtypes; i++) {
    sto->indmat[strattailattrs[i]][stratheadattrs[i]] = i;
    if(!DIRECTED && !BIPARTITE) {
      sto->indmat[stratheadattrs[i]][strattailattrs[i]] = i;
    }
  }

  // amat stores toggleability status of blocks mixing types
  int **amat = R_Calloc(nblockslevels, int *);
  amat[0] = INTEGER(getListElement(MHp->R, "amat"));
  for(int i = 1; i < nblockslevels; i++) {
    amat[i] = amat[i - 1] + nblockslevels;
  }

  // put edges in unsorted edgelists to start
  UnsrtEL **els = R_Calloc(sto->strat_nmixtypes, UnsrtEL *);
  for(int i = 0; i < sto->strat_nmixtypes; i++) {
    els[i] = UnsrtELInitialize(0, NULL, NULL, FALSE);
  }
  EXEC_THROUGH_NET_EDGES(tail, head, e, {
    int index = sto->indmat[sto->strat_vattr[tail]][sto->strat_vattr[head]];
    int allowed = amat[sto->blocks_vattr[tail]][sto->blocks_vattr[head]];
    if(index >= 0 && allowed) {
      UnsrtELInsert(tail, head, els[index]);
    }
  });
  R_Free(amat);

  // then initialize hash edgelists from unsorted edgelists
  sto->hash = R_Calloc(sto->strat_nmixtypes, HashEL *);
  for(int i = 0; i < sto->strat_nmixtypes; i++)
    sto->hash[i] = UnsrtELIntoHashEL(els[i]);
  R_Free(els);

  // initialize sampling weights for strat mixing types
  sto->original_weights = R_Calloc(sto->strat_nmixtypes, double);
  if(asInteger(getListElement(MHp->R, "empirical_flag"))) {
    // use edgecounts as weights
    for(int i = 0; i < sto->strat_nmixtypes; i++) {
      if(HashELSize(sto->hash[i]) > 0) {
        sto->original_weights[i] = HashELSize(sto->hash[i]);
      }
    }
  } else {
    // use user-supplied weights
    memcpy(sto->original_weights, 
           REAL(getListElement(MHp->R, "probvec")),
           sto->strat_nmixtypes*sizeof(double));
  }

  // determine what strat mixing types are initially toggleable
  double *currentprobvec = R_Calloc(sto->strat_nmixtypes, double);
  for(int i = 0; i < sto->strat_nmixtypes; i++) {
    // if any edges or dyads of this type are toggleable, then the mixing type is toggleable
    if(HashELSize(sto->hash[i]) > 0 || BDStratBlocksDyadCountPositive(sto->blocks, i)) {
      currentprobvec[i] = sto->original_weights[i];
      sto->current_total_weight += sto->original_weights[i];
    }
  }

  // initialize weighted sampling data structure
  sto->wtp = WtPopInitialize(sto->strat_nmixtypes,
                             currentprobvec,
                             asInteger(getListElement(MHp->R, "dyad_indep")) ? 'W' : 'B');
  R_Free(currentprobvec);

  // check degree bounds
  for(Vertex vertex = 1; vertex <= N_NODES; vertex++) {
    for(int i = 0; i < sto->bd_nlevels; i++) {
      if(DIRECTED ? (sto->indegree[i][vertex] > sto->maxin[i][vertex]
                     || sto->outdegree[i][vertex] > sto->maxout[i][vertex])
                  : (sto->indegree[i][vertex] + sto->outdegree[i][vertex] > sto->maxout[i][vertex])) {
        error("degree bound is violated by initial network; proposal cannot proceed");
      }
    }
  }
}

MH_P_FN(MH_BDStratTNT) {
  GET_STORAGE(BDStratTNTStorage, sto);

  // sample a toggleable strat mixing type on which to make a proposal
  sto->stratmixingtype = WtPopGetRand(sto->wtp);

  // number of (toggleable) edges and (toggleable) submaximal dyads of this mixing type
  int nedgestype = HashELSize(sto->hash[sto->stratmixingtype]);
  Dyad ndyadstype = BDStratBlocksDyadCount(sto->blocks, sto->stratmixingtype);

  int edgestate;
  if((unif_rand() < 0.5 && nedgestype > 0) || ndyadstype == 0) {
    // propose toggling off a random (toggleable) edge of the sampled strat mixing type
    HashELGetRand(Mtail, Mhead, sto->hash[sto->stratmixingtype]);
    edgestate = TRUE;
  } else {
    // proposed toggling a random (toggleable) submaximal dyad of the sampled strat mixing type
    BDStratBlocksGetRandWithCount(Mtail, Mhead, sto->blocks, sto->stratmixingtype, ndyadstype);
    edgestate = IS_OUTEDGE(*Mtail, *Mhead);
  }

  // determine if tail and/or head will change degree maximality status
  // if the proposed toggle is accepted
  int tailattr = sto->bd_vattr[*Mtail];
  int headattr = sto->bd_vattr[*Mhead];
  sto->tailmaxl = (DIRECTED ? sto->outdegree[headattr][*Mtail]
                            : sto->indegree[headattr][*Mtail]
                              + sto->outdegree[headattr][*Mtail])
                  == sto->maxout[headattr][*Mtail] - 1 + edgestate;
  sto->headmaxl = (DIRECTED ? sto->indegree[tailattr][*Mhead]
                            : sto->indegree[tailattr][*Mhead]
                              + sto->outdegree[tailattr][*Mhead])
                  == sto->maxin[tailattr][*Mhead] - 1 + edgestate;

  // compute proposed dyad count for current mixing type (only)
  Dyad proposedndyadstype = BDStratBlocksDyadCountOnToggle(*Mtail, *Mhead, sto->blocks,
                                                           sto->stratmixingtype,
                                                           sto->tailmaxl,
                                                           sto->headmaxl);

  // determine if any other strat mixing types will have their toggleability
  // status change if the proposed toggle is accepted
  ComputeChangesToToggleability(Mtail, Mhead, sto);

  // the rationale for the forward proposal probability (conditional on having
  // sampled the current strat mixing type) is as follows:
  // - if the sampled dyad is an edge, then we could have sampled it through
  //   the "random edge" branch above, which we enter with probability 1 if there
  //   are no toggleable, submaximal dyads of the current strat mixing type, and
  //   with probability 1/2 if there are toggleable, submaximal dyads of the current
  //   strat mixing type; additionally, if there are toggleable, submaximal dyads
  //   of the current strat mixing type *and* both tail and head are submaximal
  //   in the current network, then we could also have selected the sampled dyad
  //   through the "random dyad" branch above, which we enter with probability
  //   1/2 (since the sampled dyad is an edge and thus there are (toggleable) edges
  //   of the current strat mixing type)
  // - if the sampled dyad is not an edge, then we could only have sampled it
  //   through the "random dyad" branch, which we enter with probability 1 if
  //   there are no (toggleable) edges of the current strat mixing type, and with
  //   probability 1/2 if there are edges of the current strat mixing type
  double forward = edgestate
                   ? (ndyadstype == 0
                      ? 1.0/nedgestype
                      : 0.5/nedgestype + (sto->tailmaxl || sto->headmaxl
                                          ? 0.0
                                          : 0.5/ndyadstype))
                   : (nedgestype == 0 ? 1.0 : 0.5)/ndyadstype;

  // the backward proposal probability is basically an inverted (with respect to
  // edgestate) form of the forward proposal probability, with edge and dyad counts
  // updated to reflect the state of the network that will result if the current
  // proposal is accepted
  double backward = edgestate
                    ? (nedgestype == 1 ? 1.0 : 0.5)/proposedndyadstype
                    : (proposedndyadstype == 0
                       ? 1.0/(nedgestype + 1)
                       : 0.5/(nedgestype + 1) + (sto->tailmaxl || sto->headmaxl
                                                 ? 0.0
                                                 : 0.5/proposedndyadstype));

  // the probability to select the current strat mixing type (which is necessarily
  // toggleable in both the current and proposed networks) is inversely proportional
  // to the total weight of toggleable strat mixing types; this total weight can differ
  // between the current and proposed networks, which we account for here
  double prob_weight = sto->current_total_weight/sto->proposed_total_weight;

  MHp->logratio = log(prob_weight*backward/forward);
}

MH_U_FN(Mu_BDStratTNT) {
  GET_STORAGE(BDStratTNTStorage, sto);

  int tailattr = sto->bd_vattr[tail];
  int headattr = sto->bd_vattr[head];

  if(sto->CD) {
    // if we are running the contrastive divergence algorithm, we can't assume
    // that the tail and head passed to this U_FN have just come from the P_FN,
    // so we recompute some things here that are typically handled in the P_FN
    sto->stratmixingtype = sto->indmat[sto->strat_vattr[tail]][sto->strat_vattr[head]];

    sto->tailmaxl = (DIRECTED ? sto->outdegree[headattr][tail]
                              : sto->indegree[headattr][tail] + sto->outdegree[headattr][tail])
                     == sto->maxout[headattr][tail] - 1 + edgestate;
    sto->headmaxl = (DIRECTED ? sto->indegree[tailattr][head]
                              : sto->indegree[tailattr][head] + sto->outdegree[tailattr][head])
                     == sto->maxin[tailattr][head] - 1 + edgestate;

    ComputeChangesToToggleability(&tail, &head, sto);
  }

  // update degree information
  sto->indegree[tailattr][head] += edgestate ? -1 : 1;
  sto->outdegree[headattr][tail] += edgestate ? -1 : 1;

  // update edgelist
  HashELToggleKnown(tail, head, sto->hash[sto->stratmixingtype], edgestate);

  // update nodelists as needed
  BDNodeListsToggleIf(tail, head, sto->lists, sto->tailmaxl, sto->headmaxl);

  // if any strat mixing types have changed toggleability status,
  // update prob info accordingly
  if(sto->strat_nmixtypestoupdate > 0) {
    sto->current_total_weight = sto->proposed_total_weight;
    for(int i = 0; i < sto->strat_nmixtypestoupdate; i++) {
      WtPopSetWt(sto->strat_mixtypestoupdate[i],
                 edgestate ? sto->original_weights[sto->strat_mixtypestoupdate[i]] : 0,
                 sto->wtp);
    }
  }
}

MH_F_FN(Mf_BDStratTNT) {
  // free memory allocated in I_FN
  GET_STORAGE(BDStratTNTStorage, sto);

  BDNodeListsDestroy(sto->lists);
  BDStratBlocksDestroy(sto->blocks);

  for(int i = 0; i < sto->strat_nmixtypes; i++) {
    HashELDestroy(sto->hash[i]);
  }
  R_Free(sto->hash);

  for(int i = 0; i < sto->strat_nlevels; i++) {
    R_Free(sto->indmat[i]);
  }
  R_Free(sto->indmat);

  R_Free(sto->original_weights);

  R_Free(sto->strat_mixtypestoupdate);

  WtPopDestroy(sto->wtp);

  R_Free(sto->maxout);
  if(DIRECTED) {
    R_Free(sto->maxin);
  }

  for(int i = 0; i < sto->bd_nlevels; i++) {
    R_Free(sto->indegree[i]);
    R_Free(sto->outdegree[i]);
  }
  R_Free(sto->indegree);
  R_Free(sto->outdegree);
  // MHp->storage itself should be R_Freed by MHProposalDestroy
}

/********************
   void MH_TNT10
   Attempts to do 10 TNT steps at once, but this seems flawed currently
   because it does not correctly update network quantities like nedges
   after each of the 10 proposed toggles.
***********************/

MH_I_FN(Mi_TNT10){
  MH_STORAGE = DegreeBoundInitializeR(MHp->R, nwp);
  MHp->ntoggles=10;
}

MH_P_FN(MH_TNT10){
  Edge nedges=EDGECOUNT(nwp);
  const double P=0.5;
  double Q = 1-P;
  double DP = P*DYADCOUNT(nwp);
  double DO = DP/Q;


  double logratio = 0;
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

  CHECK_BD(MH_STORAGE);
  MHp->logratio += logratio;
}

MH_F_FN(Mf_TNT10){
  DegreeBoundDestroy(MH_STORAGE);
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
MH_I_FN(Mi_ConstantEdges){
  INIT_DYADGEN_AND_DEGREE_BOUND_AND_MODEL(TRUE);
  Edge nedges = DyadGenEdgecount(storage->gen);
  MHp->ntoggles = nedges>0 && nedges<storage->gen->ndyads && storage->gen->ndyads>=2 ? 2 : MH_FAILED;
}

MH_P_FN(MH_ConstantEdges){  
  GET_STORAGE(StoreDyadGenAndDegreeBoundAndModel, storage);

  /* First, select edge at random */
  DyadGenRandEdge(Mtail, Mhead, storage->gen);
  /* Second, select non-edge at random */
  DyadGenRandNonedge(Mtail+1, Mhead+1, storage->gen);

  CHECK_BD(storage->bd);
  CHECK_CHANGESTATS(storage->m);
}

MH_F_FN(Mf_ConstantEdges){
  DESTROY_DYADGEN_AND_DEGREE_BOUND_AND_MODEL;
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
MH_I_FN(Mi_CondOutDegreeDist){
  MH_STORAGE = DegreeBoundInitializeR(MHp->R, nwp);
  MHp->ntoggles=2;
}

MH_P_FN(MH_CondOutDegreeDist){  
  int noutedge=0, k, fvalid=0;
  int k0, k1;
  int trynode;
  Vertex e, alter, tail=0, head, head1;
  
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
  
  if(trynode==1500 || !CheckTogglesValid(MH_STORAGE, MHp, nwp)){
      Mtail[0] = 1;
      Mhead[0] = 2;
      Mtail[1] = 1;
      Mhead[1] = 2;
  }
  

}

MH_F_FN(Mf_CondOutDegreeDist){
  DegreeBoundDestroy(MH_STORAGE);
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


MH_I_FN(Mi_NodePairedTiesToggles){
  MH_STORAGE = DegreeBoundInitializeR(MHp->R, nwp);
  MHp->ntoggles = N_NODES;
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
  if (!CheckTogglesValid(MH_STORAGE, MHp, nwp))
    {
      *Mtail = *Mhead = 0;
    }
}

MH_F_FN(Mf_NodePairedTiesToggles){
  DegreeBoundDestroy(MH_STORAGE);
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

  edges = (Vertex *) R_Calloc(N_NODES+1, Vertex);
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
  R_Free(edges);
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
  edges1 = (Vertex *) R_Calloc(N_NODES+1, Vertex);
  edges2 = (Vertex *) R_Calloc(N_NODES+1, Vertex);
  
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
  R_Free(edges1);
  R_Free(edges2);
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
  outedges = (Vertex *) R_Calloc(N_NODES+1, Vertex);
  inedges = (Vertex *) R_Calloc(N_NODES+1, Vertex);
  
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
  
  R_Free(outedges);
  R_Free(inedges);
  
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


MH_I_FN(Mi_ConstrainedNodePairedTiesToggles){
  MH_STORAGE = DegreeBoundInitializeR(MHp->R, nwp);
  MHp->ntoggles=N_NODES;
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
  if (!CheckConstrainedTogglesValid(MH_STORAGE, MHp, nwp))
    {
      *Mtail = *Mhead = 0;
    }
}

MH_F_FN(Mf_ConstrainedNodePairedTiesToggles){
  DegreeBoundDestroy(MH_STORAGE);
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

  edges = (Vertex *) R_Calloc(N_NODES+1, Vertex);
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
  R_Free(edges);
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
MH_I_FN(Mi_ConstrainedTwoRandomToggles){
  MH_STORAGE = DegreeBoundInitializeR(MHp->R, nwp);
  MHp->ntoggles=2;
}

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
  
  if (!CheckConstrainedTogglesValid(MH_STORAGE, MHp, nwp))
    {
      Mtail[0] = Mhead[0] = 0;
      Mtail[1] = Mhead[1] = 0;
    }  
}

MH_F_FN(Mf_ConstrainedTwoRandomToggles){
  DegreeBoundDestroy(MH_STORAGE);
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
  edges1 = (Vertex *) R_Calloc(N_NODES+1, Vertex);
  edges2 = (Vertex *) R_Calloc(N_NODES+1, Vertex);
  
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
    R_Free(edges1);
    R_Free(edges2);
      }
  if (tail2 > head2)
    {
      Mtail[1] = head2;
      Mhead[1] = tail2;
    }else{
      Mtail[1] = tail2;
      Mhead[1] = head2;
    }
  R_Free(edges1);
  R_Free(edges2);
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

  edges1 = (Vertex *) R_Calloc(N_NODES+1, Vertex);
  edges2 = (Vertex *) R_Calloc(N_NODES+1, Vertex);

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
  R_Free(edges1);
  R_Free(edges2);
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


