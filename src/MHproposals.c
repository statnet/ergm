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
#include "ergm_nodelist_dyad_sampler.h"
#include "ergm_BDStrat_proposals.h"

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

  sto->bound = asInteger(getListElement(MHp->R, "bound"));
  sto->nmixtypes = asInteger(getListElement(MHp->R, "nmixtypes"));
  sto->nstratlevels = asInteger(getListElement(MHp->R, "nattrcodes"));
  sto->nbdlevels = asInteger(getListElement(MHp->R, "bd_levels"));
  sto->strattailtypes = INTEGER(getListElement(MHp->R, "strattailattrs"));
  sto->stratheadtypes = INTEGER(getListElement(MHp->R, "stratheadattrs"));

  sto->mixtypestoupdate = Calloc(sto->nmixtypes, int);
  
  sto->strat_vattr = INTEGER(getListElement(MHp->R, "strat_vattr"));
  sto->strat_vattr--; // so node indices line up correctly

  sto->bd_vattr = INTEGER(getListElement(MHp->R, "bd_vattr"));
  sto->bd_vattr--; // so node indices line up correctly  
  
  
  sto->BDtypesbyStrattype = INTEGER(getListElement(MHp->R, "BDtypesbyStrattype"));
  
  sto->BDtailsbyStrattype = Calloc(sto->nmixtypes, int *);
  sto->BDheadsbyStrattype = Calloc(sto->nmixtypes, int *);
  
  sto->BDtailsbyStrattype[0] = INTEGER(getListElement(MHp->R, "BDtailsbyStrattype"));
  sto->BDheadsbyStrattype[0] = INTEGER(getListElement(MHp->R, "BDheadsbyStrattype"));
  
  for(int i = 1; i < sto->nmixtypes; i++) {
    sto->BDtailsbyStrattype[i] = sto->BDtailsbyStrattype[i - 1] + sto->BDtypesbyStrattype[i - 1];
    sto->BDheadsbyStrattype[i] = sto->BDheadsbyStrattype[i - 1] + sto->BDtypesbyStrattype[i - 1];
  }
  
  
  sto->els = Calloc(sto->nmixtypes, UnsrtEL *);
  for(int i = 0; i < sto->nmixtypes; i++) {
    sto->els[i] = UnsrtELInitialize(0, NULL, NULL, FALSE);
  }
      
  sto->indmat = Calloc(sto->nstratlevels, int *);
  sto->indmat[0] = INTEGER(getListElement(MHp->R, "indmat"));
  for(int i = 1; i < sto->nstratlevels; i++) {
    sto->indmat[i] = sto->indmat[i - 1] + sto->nstratlevels;
  }
  
  sto->currentsubmaxledgestype = Calloc(sto->nmixtypes, int);  
  EXEC_THROUGH_NET_EDGES(tail, head, e, {
    int index = sto->indmat[sto->strat_vattr[tail]][sto->strat_vattr[head]];
    if(index >= 0) {
      UnsrtELInsert(tail, head, sto->els[index]);
      if(IN_DEG[tail] + OUT_DEG[tail] < sto->bound && IN_DEG[head] + OUT_DEG[head] < sto->bound) {
        sto->currentsubmaxledgestype[index]++;
      }
    }  
  });  
    
  int *nodecounts = INTEGER(getListElement(MHp->R, "nodecountsbypairedcode"));
  
  sto->nodesvec = Calloc(sto->nstratlevels, Vertex **);
  for(int i = 0; i < sto->nstratlevels; i++) {
    sto->nodesvec[i] = Calloc(sto->nbdlevels, Vertex *);
    for(int j = 0; j < sto->nbdlevels; j++) {
      sto->nodesvec[i][j] = Calloc(nodecounts[i*sto->nbdlevels + j], Vertex);
    }
  }
  
  sto->attrcounts = Calloc(sto->nstratlevels, int *);
  for(int i = 0; i < sto->nstratlevels; i++) {
    sto->attrcounts[i] = (int *)Calloc(sto->nbdlevels, int);
  }

  sto->nodepos = Calloc(N_NODES + 1, int);  

  for(Vertex vertex = 1; vertex <= N_NODES; vertex++) {
    if(IN_DEG[vertex] + OUT_DEG[vertex] < sto->bound) {
      // add vertex to the submaximal list corresponding to its attribute type
      sto->nodesvec[sto->strat_vattr[vertex]][sto->bd_vattr[vertex]][sto->attrcounts[sto->strat_vattr[vertex]][sto->bd_vattr[vertex]]] = vertex;
      sto->nodepos[vertex] = sto->attrcounts[sto->strat_vattr[vertex]][sto->bd_vattr[vertex]];      
      sto->attrcounts[sto->strat_vattr[vertex]][sto->bd_vattr[vertex]]++;
    }
  }
  
  sto->originalprobvec = Calloc(sto->nmixtypes, double);
  int empirical_flag = asInteger(getListElement(MHp->R, "empirical_flag"));
  if(empirical_flag) {
    for(int i = 0; i < sto->nmixtypes; i++) {
      if(sto->els[i]->nedges > 0) {
        sto->originalprobvec[i] = sto->els[i]->nedges;
        sto->nmixtypes_max++;
      }
    }
  } else {
    memcpy(sto->originalprobvec, REAL(getListElement(MHp->R, "probvec")), sto->nmixtypes*sizeof(double));
    sto->nmixtypes_max = sto->nmixtypes;
  }

  // determine what mixing types are initially toggleable 
  double *currentprobvec = Calloc(sto->nmixtypes, double);  
  for(int i = 0; i < sto->nmixtypes; i++) {
    // if any edges or dyads of this type are toggleable
    if(sto->els[i]->nedges > 0 || NodeListDyadCountPositive(sto->attrcounts[sto->strattailtypes[i]], sto->attrcounts[sto->stratheadtypes[i]], sto->BDtailsbyStrattype[i], sto->BDheadsbyStrattype[i], sto->BDtypesbyStrattype[i], sto->strattailtypes[i] == sto->stratheadtypes[i])) {
      currentprobvec[i] = sto->originalprobvec[i];
      sto->currentcumprob += sto->originalprobvec[i];
      sto->nmixtypes_toggleable++;
    }
  }
    
  sto->wtp = WtPopInitialize(sto->nmixtypes, currentprobvec);
  Free(currentprobvec);

  // zero proposal probability is an error
  if(WtPopSumWts(sto->wtp) == 0) {
    MHp->ntoggles = MH_FAILED;
    return;
  }
}

MH_P_FN(MH_BDStratTNT) {
  GET_STORAGE(BDStratTNTStorage, sto);

  // sample a toggleable strat mixing type on which to make a proposal
  int strat_i = WtPopGetRand(sto->wtp);

  // record the mixing type of the toggle, in case it's needed in the U function later
  sto->stratmixingtype = strat_i;    

  int strattailtype = sto->strattailtypes[strat_i];
  int stratheadtype = sto->stratheadtypes[strat_i];
    
  // number of edges of this mixing type
  int nedgestype = sto->els[strat_i]->nedges;
  
  Dyad ndyadstype = NodeListDyadCount(sto->attrcounts[strattailtype], sto->attrcounts[stratheadtype], sto->BDtailsbyStrattype[strat_i], sto->BDheadsbyStrattype[strat_i], sto->BDtypesbyStrattype[strat_i], strattailtype == stratheadtype, DIRECTED);
  
  int edgeflag;
  
  if((unif_rand() < 0.5 && nedgestype > 0) || ndyadstype == 0) {
    // propose toggling off an existing edge of strat mixing type strat_i
    UnsrtELGetRand(Mtail, Mhead, sto->els[strat_i]);    
    edgeflag = TRUE;
  } else {
    // select a random BD toggleable dyad of strat mixing type strat_i and propose toggling it
    GetRandDyadFromLists(Mtail, // tail
                         Mhead, // head
                         sto->nodesvec[strattailtype], // tails
                         sto->nodesvec[stratheadtype], // heads
                         sto->BDtailsbyStrattype[strat_i], // tailattrs
                         sto->BDheadsbyStrattype[strat_i], // headattrs
                         sto->attrcounts[strattailtype], // tailcounts
                         sto->attrcounts[stratheadtype], // headcounts
                         sto->BDtypesbyStrattype[strat_i], // length
                         ndyadstype, // dyadcount
                         strattailtype == stratheadtype, // diagonal
                         DIRECTED); // directed; always FALSE in BDStratTNT

    // now check if the dyad we drew is already an edge or not
    if(IS_OUTEDGE(Mtail[0],Mhead[0])) {
      // must resample to know edge index; we will fix strat and bd types below;
      // indices within submaximal lists will not be used in the U_FN since this is an off-toggle,
      // so if they're wrong, that's fine
      UnsrtELGetRand(Mtail, Mhead, sto->els[strat_i]);
      edgeflag = TRUE;
    } else {
      edgeflag = FALSE;
    }
  }

  sto->strattailtype = sto->strat_vattr[Mtail[0]];
  sto->stratheadtype = sto->strat_vattr[Mhead[0]];
  
  sto->bdtailtype = sto->bd_vattr[Mtail[0]];
  sto->bdheadtype = sto->bd_vattr[Mhead[0]];
  
  sto->tailmaxl = IN_DEG[Mtail[0]] + OUT_DEG[Mtail[0]] == sto->bound - 1 + edgeflag;
  sto->headmaxl = IN_DEG[Mhead[0]] + OUT_DEG[Mhead[0]] == sto->bound - 1 + edgeflag;
    
    
  // temporarily set tail and head toggleability to what it would be in the proposed network
  if(sto->tailmaxl) {
    NodeListToggleKnown(*Mtail, sto->nodesvec[sto->strattailtype][sto->bdtailtype], sto->nodepos, sto->attrcounts[sto->strattailtype] + sto->bdtailtype, !edgeflag);
  }
  if(sto->headmaxl) {
    NodeListToggleKnown(*Mhead, sto->nodesvec[sto->stratheadtype][sto->bdheadtype], sto->nodepos, sto->attrcounts[sto->stratheadtype] + sto->bdheadtype, !edgeflag);
  }

  // compute proposed dyad count for current mixing type (only)
  Dyad proposeddyadstype = NodeListDyadCount(sto->attrcounts[strattailtype], sto->attrcounts[stratheadtype], sto->BDtailsbyStrattype[strat_i], sto->BDheadsbyStrattype[strat_i], sto->BDtypesbyStrattype[strat_i], strattailtype == stratheadtype, DIRECTED);
    
  // here we compute the proposedcumprob, checking only those
  // mixing types that can be influenced by toggles made on 
  // the current mixing type
  sto->proposedcumprob = sto->currentcumprob;
  sto->nmixtypestoupdate = 0; // reset counter
  // avoid these somewhat expensive checks in the typical case
  // where you have enough submaximal nodes that you cannot
  // be exhausting any mixing types of toggleable dyads
  if(sto->attrcounts[sto->strattailtype][sto->bdtailtype] <= 2 || sto->attrcounts[sto->stratheadtype][sto->bdheadtype] <= 2) {
    
    // how many strat types do we need to check?
    int ntocheck = (sto->strattailtype == sto->stratheadtype) ? sto->nstratlevels : 2*sto->nstratlevels;

    for(int i = 0; i < ntocheck; i++) {
      // find the index of the i'th strat type we need to check, by looking it up in the indmat
      int infl_i = sto->indmat[i < sto->nstratlevels ? sto->strattailtype : i - sto->nstratlevels][i < sto->nstratlevels ? i : sto->stratheadtype];

      // if this strat type is not included in the proposal, or is the same as the strat type of the proposed toggle,
      // then it cannot change toggleability status, so skip it
      if(infl_i < 0 || infl_i == strat_i) {
        continue;
      }
      
      // can we toggle this mixing type in the current network?
      int toggle_curr = WtPopGetWt(infl_i, sto->wtp) > 0;
      
      // will we be able to toggle this mixing type in the proposed network? 
      int toggle_prop = sto->els[infl_i]->nedges > 0 || NodeListDyadCountPositive(sto->attrcounts[sto->strattailtypes[infl_i]], sto->attrcounts[sto->stratheadtypes[infl_i]], sto->BDtailsbyStrattype[infl_i], sto->BDheadsbyStrattype[infl_i], sto->BDtypesbyStrattype[infl_i], sto->strattailtypes[infl_i] == sto->stratheadtypes[infl_i]);
      
      // will there be a change in toggleability status?
      int change = toggle_curr - toggle_prop;

      // if so, take this into account      
      if(change) {
        sto->proposedcumprob -= change*sto->originalprobvec[infl_i];
        sto->mixtypestoupdate[sto->nmixtypestoupdate] = infl_i;
        sto->nmixtypestoupdate++;        
      }
    }
  }
  
  // restore tail and head toggleability to their current status
  if(sto->tailmaxl) {
    NodeListToggleKnown(*Mtail, sto->nodesvec[sto->strattailtype][sto->bdtailtype], sto->nodepos, sto->attrcounts[sto->strattailtype] + sto->bdtailtype, edgeflag);
  }
  if(sto->headmaxl) {
    NodeListToggleKnown(*Mhead, sto->nodesvec[sto->stratheadtype][sto->bdheadtype], sto->nodepos, sto->attrcounts[sto->stratheadtype] + sto->bdheadtype, edgeflag);
  }


  int delta = edgeflag ? +1 : -1;
  // for the logratio, we will need to know the number of submaxl edges in the proposed network
  int proposedsubmaxledgestype = sto->currentsubmaxledgestype[strat_i];
  
  // if we are adding an edge that will be submaximal in the post-toggle 
  // network, then increment proposedsubmaxledgestype for this particular edge
  if(!edgeflag && !sto->tailmaxl && !sto->headmaxl) {
    proposedsubmaxledgestype++;
  }
  
  // if we are removing an edge that is submaximal in the current
  // network, decrement proposedsubmaxledgestype for this particular edge
  if(edgeflag && !sto->tailmaxl && !sto->headmaxl) {
    proposedsubmaxledgestype--;
  }

  Edge e;
  Vertex v;

  // if tail will change maximality on toggle, then adjust
  // proposedsubmaxledgestype for all edges between tail and
  // a submaximal neighbor v with the edge between tail and v
  // having the mixing type strat_i, taking care not to count head,
  // since that was handled separately above
  if(sto->tailmaxl) {
    STEP_THROUGH_OUTEDGES(Mtail[0], e, v) {
      if(v != Mhead[0] && IN_DEG[v] + OUT_DEG[v] < sto->bound && sto->indmat[sto->strattailtype][sto->strat_vattr[v]] == strat_i) {
        proposedsubmaxledgestype += delta;
      }
    }
    STEP_THROUGH_INEDGES(Mtail[0], e, v) {
      if(IN_DEG[v] + OUT_DEG[v] < sto->bound && sto->indmat[sto->strattailtype][sto->strat_vattr[v]] == strat_i) {
        proposedsubmaxledgestype += delta;
      }
    }
  }
  
  // ditto head
  if(sto->headmaxl) {
    STEP_THROUGH_OUTEDGES(Mhead[0], e, v) {
      if(IN_DEG[v] + OUT_DEG[v] < sto->bound && sto->indmat[sto->stratheadtype][sto->strat_vattr[v]] == strat_i) {
        proposedsubmaxledgestype += delta;
      }
    }
    STEP_THROUGH_INEDGES(Mhead[0], e, v) {
      if(v != Mtail[0] && IN_DEG[v] + OUT_DEG[v] < sto->bound && sto->indmat[sto->stratheadtype][sto->strat_vattr[v]] == strat_i) {
        proposedsubmaxledgestype += delta;
      }
    }
  }
  
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
    MHp->logratio = log(prob_weight*(((nedgestype == 1 ? 1.0 : 0.5)/proposeddyadstype))/(((ndyadstype == 0 ? 1.0/nedgestype : (0.5/nedgestype)*(1 + (double)sto->currentsubmaxledgestype[strat_i]/ndyadstype)))));
  } else {
    MHp->logratio = log(prob_weight*((proposeddyadstype == 0 ? 1.0/(nedgestype + 1) : (0.5/(nedgestype + 1))*(1 + (double)proposedsubmaxledgestype/proposeddyadstype))/((nedgestype == 0 ? 1.0 : 0.5)/ndyadstype)));
  }
}

MH_U_FN(Mu_BDStratTNT) {
  GET_STORAGE(BDStratTNTStorage, sto);

  // update edgelist
  UnsrtELToggleKnown(tail, head, sto->els[sto->stratmixingtype], edgeflag);

  // update nodelists
  if(sto->tailmaxl) {
    NodeListToggleKnown(tail, sto->nodesvec[sto->strattailtype][sto->bdtailtype], sto->nodepos, sto->attrcounts[sto->strattailtype] + sto->bdtailtype, !edgeflag);
  }
  if(sto->headmaxl) {
    NodeListToggleKnown(head, sto->nodesvec[sto->stratheadtype][sto->bdheadtype], sto->nodepos, sto->attrcounts[sto->stratheadtype] + sto->bdheadtype, !edgeflag);
  }

  // if any strat mixing types have changed toggleability status, update prob info accordingly
  if(sto->nmixtypestoupdate > 0) {
    sto->currentcumprob = sto->proposedcumprob;

    if(edgeflag) {
      sto->nmixtypes_toggleable += sto->nmixtypestoupdate;
      for(int i = 0; i < sto->nmixtypestoupdate; i++) {
        WtPopSetWt(sto->mixtypestoupdate[i], sto->originalprobvec[sto->mixtypestoupdate[i]], sto->wtp);          
      }
    } else {
      sto->nmixtypes_toggleable -= sto->nmixtypestoupdate;
      for(int i = 0; i < sto->nmixtypestoupdate; i++) {
        WtPopSetWt(sto->mixtypestoupdate[i], 0, sto->wtp);          
      }
    }
  }
  
  // if we are adding an edge that will be submaximal in the post-toggle 
  // network, then increment currentsubmaxledgestype for this particular edge
  if(!edgeflag && !sto->tailmaxl && !sto->headmaxl) {
    sto->currentsubmaxledgestype[sto->stratmixingtype]++;
  }
  
  // if we are removing an edge that is submaximal in the current
  // network, decrement currentsubmaxledgestype for this particular edge
  if(edgeflag && !sto->tailmaxl && !sto->headmaxl) {
    sto->currentsubmaxledgestype[sto->stratmixingtype]--;
  }

  int delta = edgeflag ? +1 : -1;

  Edge e;
  Vertex v;

  // if tail will change maximality on toggle, then adjust
  // currentsubmaxledgestype for all edges between tail and
  // a submaximal neighbor v, taking care not to count head,
  // since that was handled separately above
  if(sto->tailmaxl) {
    STEP_THROUGH_OUTEDGES(tail, e, v) {
      if(v != head && IN_DEG[v] + OUT_DEG[v] < sto->bound && sto->indmat[sto->strattailtype][sto->strat_vattr[v]] >= 0) {
        sto->currentsubmaxledgestype[sto->indmat[sto->strattailtype][sto->strat_vattr[v]]] += delta;
      }
    }
    STEP_THROUGH_INEDGES(tail, e, v) {
      if(IN_DEG[v] + OUT_DEG[v] < sto->bound && sto->indmat[sto->strattailtype][sto->strat_vattr[v]] >= 0) {
        sto->currentsubmaxledgestype[sto->indmat[sto->strattailtype][sto->strat_vattr[v]]] += delta;
      }
    }
  }
  
  // ditto head
  if(sto->headmaxl) {
    STEP_THROUGH_OUTEDGES(head, e, v) {
      if(IN_DEG[v] + OUT_DEG[v] < sto->bound && sto->indmat[sto->stratheadtype][sto->strat_vattr[v]] >= 0) {
        sto->currentsubmaxledgestype[sto->indmat[sto->stratheadtype][sto->strat_vattr[v]]] += delta;
      }
    }
    STEP_THROUGH_INEDGES(head, e, v) {
      if(v != tail && IN_DEG[v] + OUT_DEG[v] < sto->bound && sto->indmat[sto->stratheadtype][sto->strat_vattr[v]] >= 0) {
        sto->currentsubmaxledgestype[sto->indmat[sto->stratheadtype][sto->strat_vattr[v]]] += delta;
      }
    }
  }
}

MH_F_FN(Mf_BDStratTNT) {
  // Free all the things
  GET_STORAGE(BDStratTNTStorage, sto);
  
  for(int i = 0; i < sto->nstratlevels; i++) {
    Free(sto->attrcounts[i]);
  }
  Free(sto->attrcounts);

  for(int i = 0; i < sto->nstratlevels; i++) {
    for(int j = 0; j < sto->nbdlevels; j++) {
      Free(sto->nodesvec[i][j]);
    }
    Free(sto->nodesvec[i]);
  }
  Free(sto->nodesvec);
  Free(sto->nodepos);
  
  for(int i = 0; i < sto->nmixtypes; i++) {
    UnsrtELDestroy(sto->els[i]);
  }
  Free(sto->els);

  Free(sto->currentsubmaxledgestype);
  Free(sto->indmat);

  Free(sto->BDtailsbyStrattype);
  Free(sto->BDheadsbyStrattype);

  Free(sto->originalprobvec);

  Free(sto->mixtypestoupdate);  

  WtPopDestroy(sto->wtp);
  // MHp->storage itself should be Freed by MHProposalDestroy
}


/********************
    MH_BDTNT

This proposal handles a constant outdegree upper bound and uses a
whitelist of blocks (sets of dyads defined by combinations of vertex
attributes) so that specific combinations can be forbidden.
********************/

// struct definition in ergm_BDStrat_proposals.h

MH_I_FN(Mi_BDTNT) {
  // process the inputs and initialize all the node lists in storage; set MHp->ntoggles to 1
  MHp->ntoggles = 1;

  // used during initialization only; not retained in storage
  int *nodecountsbycode = INTEGER(getListElement(MHp->R, "nodecountsbycode")); // Number of nodes of each type.

  ALLOC_STORAGE(1, BDTNTStorage, sto);
  sto->bound = asInteger(getListElement(MHp->R, "bound")); // As in struct.
  sto->nlevels = asInteger(getListElement(MHp->R, "nlevels")); // Number of distinct types of types of vertex.
  sto->nmixtypes = asInteger(getListElement(MHp->R, "nmixtypes")); // As in struct.
  sto->tailtypes = INTEGER(getListElement(MHp->R, "allowed.tails")); // As in struct.
  sto->headtypes = INTEGER(getListElement(MHp->R, "allowed.heads")); // As in struct.
  
  sto->vattr = INTEGER(getListElement(MHp->R, "nodecov")); // As in struct.
  sto->vattr--; // so node indices line up correctly
  
  sto->nodesvec = Calloc(sto->nlevels, Vertex *); // As in struct.
  sto->nodepos = Calloc(N_NODES + 1, int); // As in struct.
  sto->attrcounts = Calloc(sto->nlevels, int); // As in struct.
  // make room for maximum number of nodes of each type
  for(int i = 0; i < sto->nlevels; i++) {
    sto->nodesvec[i] = Calloc(nodecountsbycode[i], Vertex);
  }

  // Populate the list of submaximal vertices by checking whether a given vertex has a maximal degree.
  for(Vertex vertex = 1; vertex <= N_NODES; vertex++) {
    if(IN_DEG[vertex] + OUT_DEG[vertex] < sto->bound) {
      // add vertex to the submaximal list corresponding to its attribute type
      sto->nodesvec[sto->vattr[vertex]][sto->attrcounts[sto->vattr[vertex]]] = vertex;
      sto->nodepos[vertex] = sto->attrcounts[sto->vattr[vertex]];
      sto->attrcounts[sto->vattr[vertex]]++;
    }
  }
  
  // Construct and populate the list of edges. (May be obviated by more efficient network sampling.)
  sto->edgelist = UnsrtELInitialize(0, NULL, NULL, FALSE);
  EXEC_THROUGH_NET_EDGES(tail, head, e, {
    UnsrtELInsert(tail, head, sto->edgelist);
    if(IN_DEG[tail] + OUT_DEG[tail] < sto->bound && IN_DEG[head] + OUT_DEG[head] < sto->bound) {
      sto->currentsubmaxledges++;
    }
  });
  
  // count number of "BD-toggleable" dyads in current network
  sto->currentdyads = NodeListDyadCount(sto->attrcounts, sto->attrcounts, sto->tailtypes, sto->headtypes, sto->nmixtypes, TRUE, DIRECTED);

  // if we cannot toggle any edges or dyads, error
  if(EDGECOUNT(nwp) == 0 && sto->currentdyads == 0) {
    MHp->ntoggles = MH_FAILED;
    return;
  }  
}

MH_P_FN(MH_BDTNT) {    
  GET_STORAGE(BDTNTStorage, sto);

  int nedges = EDGECOUNT(nwp);
  
  int edgeflag;

  // if currentdyads == 0, we *must* propose toggling off an existing edge;
  // the case nedges == 0 && currentdyads == 0 was excluded during initialization,
  // and we cannot end up in that case if we don't start in that case
  // (assuming the initial network is valid)  
  if((unif_rand() < 0.5 && nedges > 0) || (sto->currentdyads == 0)) {
    // select an existing edge at random, and propose toggling it off
    UnsrtELGetRand(Mtail, Mhead, sto->edgelist);
        
    edgeflag = TRUE;
  } else {
    // select a BD-toggleable dyad and propose toggling it

    // Note that the dyad block selector cannot select an edge
    // incident on a maximal node; but the edge "reselection" below
    // will be able to select it with equal probability to the others,
    // equalising its chances of being selected, and requiring only a
    // "marginal" adjustment to the acceptance probability.
    GetRandDyadFromLists(Mtail, // tail
                         Mhead, // head
                         sto->nodesvec, // tails
                         sto->nodesvec, // heads
                         sto->tailtypes, // tailattrs
                         sto->headtypes, // headattrs
                         sto->attrcounts, // tailcounts
                         sto->attrcounts, // headcounts
                         sto->nmixtypes, // length
                         sto->currentdyads, // dyadcount
                         TRUE, // diagonal; no higher level types here, so always TRUE
                         DIRECTED); // directed; always FALSE in BDTNT

    edgeflag = IS_OUTEDGE(Mtail[0],Mhead[0]);
    // Resample the edge to make it known to the unsorted edgelist, if necessary
    if(edgeflag) UnsrtELGetRand(Mtail, Mhead, sto->edgelist);
  }
  
  sto->tailtype = sto->vattr[Mtail[0]];
  sto->headtype = sto->vattr[Mhead[0]];    
  
  sto->tailmaxl = IN_DEG[Mtail[0]] + OUT_DEG[Mtail[0]] == sto->bound - 1 + edgeflag;
  sto->headmaxl = IN_DEG[Mhead[0]] + OUT_DEG[Mhead[0]] == sto->bound - 1 + edgeflag;   
  
  
  // temporarily make tail and head assume their submaximal status in the proposed network
  // so we can easily compute the number of "BD toggleable dyads" in the proposed network
  if(sto->tailmaxl) {
    NodeListToggleKnown(*Mtail, sto->nodesvec[sto->tailtype], sto->nodepos, sto->attrcounts + sto->tailtype, !edgeflag);
  }
  if(sto->headmaxl) {
    NodeListToggleKnown(*Mhead, sto->nodesvec[sto->headtype], sto->nodepos, sto->attrcounts + sto->headtype, !edgeflag);    
  }
  // the count of dyads that can be toggled in the "GetRandBDDyad" branch,
  // in the proposed network
  sto->proposeddyads = NodeListDyadCount(sto->attrcounts, sto->attrcounts, sto->tailtypes, sto->headtypes, sto->nmixtypes, TRUE, DIRECTED);
  // now restore tail and head to their current state, since we won't necesssarily accept this proposed toggle
  if(sto->tailmaxl) {
    NodeListToggleKnown(*Mtail, sto->nodesvec[sto->tailtype], sto->nodepos, sto->attrcounts + sto->tailtype, edgeflag);
  }
  if(sto->headmaxl) {
    NodeListToggleKnown(*Mhead, sto->nodesvec[sto->headtype], sto->nodepos, sto->attrcounts + sto->headtype, edgeflag);    
  }
  
  // calculate the number of submaximal edges in the proposed network
  int delta = edgeflag ? +1 : -1;  
  sto->proposedsubmaxledges = sto->currentsubmaxledges;
  
  // if we are adding an edge that will be submaximal in the post-toggle 
  // network, then increment proposedsubmaxledges for this particular edge
  if(!edgeflag && !sto->tailmaxl && !sto->headmaxl) {
    sto->proposedsubmaxledges++;
  }
  
  // if we are removing an edge that is submaximal in the current
  // network, decrement proposedsubmaxledges for this particular edge
  if(edgeflag && !sto->tailmaxl && !sto->headmaxl) {
    sto->proposedsubmaxledges--;
  }

  Edge e;
  Vertex v;

  // if tail will change maximality on toggle, then adjust
  // proposedsubmaxledges for all edges between tail and
  // a submaximal neighbor v, taking care not to count head,
  // since that was handled separately above
  if(sto->tailmaxl) {
    STEP_THROUGH_OUTEDGES(Mtail[0], e, v) {
      if(v != Mhead[0] && IN_DEG[v] + OUT_DEG[v] < sto->bound) {
        sto->proposedsubmaxledges += delta;
      }
    }
    STEP_THROUGH_INEDGES(Mtail[0], e, v) {
      if(IN_DEG[v] + OUT_DEG[v] < sto->bound) {
        sto->proposedsubmaxledges += delta;
      }
    }
  }
  
  // ditto head
  if(sto->headmaxl) {
    STEP_THROUGH_OUTEDGES(Mhead[0], e, v) {
      if(IN_DEG[v] + OUT_DEG[v] < sto->bound) {
        sto->proposedsubmaxledges += delta;
      }
    }
    STEP_THROUGH_INEDGES(Mhead[0], e, v) {
      if(v != Mtail[0] && IN_DEG[v] + OUT_DEG[v] < sto->bound) {
        sto->proposedsubmaxledges += delta;
      }
    }
  }
  
  // rationale for the logratio:
  //
  // there is only one way to select a non-edge, and that is via the initial dyad sample in the "GetRandBDDyad" branch above;
  // all "BD toggleable dyads" are given equal weight in the initial dyad sampling in the GetRandBDDyad branch,
  // so if one has entered the GetRandBDDyad branch, the probability to select a given non-edge is 1 over the number of 
  // BD toggleable dyads, which is sto->currentdyads in the current network and proposeddyads in the proposed network; 
  // given that the network has at least one BD toggleable dyad, the probability to enter the GetRandBDDyad branch 
  // is either 1/2 if the network has any edges or 1 if the network has no edges; this fully explains the calculation 
  // of non-edge sampling probabilities below
  //
  // there are two ways to select an edge: it can be sampled directly in the "GetRandEdge" branch above, or it can be sampled 
  // indirectly by entering the GetRandBDDyad branch, sampling a submaximal edge, and then resampling an arbitrary edge;
  // in the GetRandEdge branch, the probability to select a given edge is 1 over the number of edges in the network;
  // in the GetRandBDDyad branch, the probability to select a submaximal edge in the initial dyad sample is the number of
  // submaximal edges divided by the number of BD toggleable dyads, and given that a submaximal edge was selected as the initial
  // dyad sample, the probability to select a given edge on resampling is 1 over the number of edges; given that an edge exists
  // in the network, the probability to enter to GetRandEdge branch is either 1/2 if any BD toggleable dyads exist in the network,
  // or is 1 if no BD toggleable dyads exist in the network; the probability to enter the GetRandBDDyad branch is 1 minus the
  // probability to enter the GetRandEdge branch; this fully explains the calculation of edge sampling probabilities below
  
  if(edgeflag) {
    MHp->logratio = log(((nedges == 1 ? 1.0 : 0.5)/sto->proposeddyads)/(sto->currentdyads == 0 ? 1.0/nedges : (0.5/nedges)*(1 + ((double)sto->currentsubmaxledges/sto->currentdyads))));
  } else {
    MHp->logratio = log((sto->proposeddyads == 0 ? 1.0/(nedges + 1) : (0.5/(nedges + 1))*(1 + (double)sto->proposedsubmaxledges/sto->proposeddyads))/((nedges == 0 ? 1.0 : 0.5)/sto->currentdyads));  
  }
}

// this U_FN is called *before* the toggle is made in the network
MH_U_FN(Mu_BDTNT) {  
  GET_STORAGE(BDTNTStorage, sto);
  // update edgelist
  UnsrtELToggleKnown(tail, head, sto->edgelist, edgeflag);

  // update nodelists
  if(sto->tailmaxl) {
    NodeListToggleKnown(tail, sto->nodesvec[sto->tailtype], sto->nodepos, sto->attrcounts + sto->tailtype, !edgeflag);
  }
  if(sto->headmaxl) {
    NodeListToggleKnown(head, sto->nodesvec[sto->headtype], sto->nodepos, sto->attrcounts + sto->headtype, !edgeflag);    
  }
  
  // update current dyad count
  sto->currentdyads = sto->proposeddyads;

  // update the current submaximal edge count
  sto->currentsubmaxledges = sto->proposedsubmaxledges;
}

MH_F_FN(Mf_BDTNT) {
  // Free all the things
  GET_STORAGE(BDTNTStorage, sto);
  UnsrtELDestroy(sto->edgelist);

  Free(sto->attrcounts);
  Free(sto->nodepos);  

  for(int i = 0; i < sto->nlevels; i++) {
    Free(sto->nodesvec[i]);
  }
  Free(sto->nodesvec);
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
  int *vattr = INTEGER(getListElement(MHp->R, "nodecov"));
  vattr--; // so node indices line up correctly

  int nlevels = asInteger(getListElement(MHp->R, "nlevels"));

  ALLOC_STORAGE(1, StratTNTStorage, sto);
  sto->nmixtypes = asInteger(getListElement(MHp->R, "nmixtypes"));
  sto->tailtypes = INTEGER(getListElement(MHp->R, "tailattrs"));
  sto->headtypes = INTEGER(getListElement(MHp->R, "headattrs"));
  sto->nodecountsbycode = INTEGER(getListElement(MHp->R, "nodecountsbycode"));
  
  sto->els = Calloc(sto->nmixtypes, UnsrtEL *);
  for(int i = 0; i < sto->nmixtypes; i++) {
    sto->els[i] = UnsrtELInitialize(0, NULL, NULL, FALSE);
  }
  
  int *inputindmat = INTEGER(getListElement(MHp->R, "indmat"));    
  int **indmat = Calloc(nlevels, int *);
  indmat[0] = inputindmat;
  for(int i = 1; i < nlevels; i++) {
    indmat[i] = indmat[i - 1] + nlevels;
  }
  
  EXEC_THROUGH_NET_EDGES(tail, head, e, {
    int index = indmat[vattr[tail]][vattr[head]];
    if(index >= 0) {
      UnsrtELInsert(tail, head, sto->els[index]);
    }      
  });
  Free(indmat);
  
  sto->nodesbycode = Calloc(nlevels, Vertex *);
  sto->nodesbycode[0] = (Vertex *)INTEGER(getListElement(MHp->R, "nodeindicesbycode"));
  for(int i = 1; i < nlevels; i++) {
    sto->nodesbycode[i] = sto->nodesbycode[i - 1] + sto->nodecountsbycode[i - 1];
  }
  
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
    sto->ndyadstype[i] = NodeListDyadCount(sto->nodecountsbycode, sto->nodecountsbycode, sto->tailtypes + i, sto->headtypes + i, 1, TRUE, DIRECTED);
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
  int i = WtPopGetRand(sto->wtp);
  
  // record the mixing type of the toggle, in case it's needed in the U function later
  sto->currentmixingtype = i;    
      
  // number of edges of this mixing type
  int nedgestype = sto->els[i]->nedges;

  // number of dyads of this mixing type
  Dyad ndyadstype = sto->ndyadstype[i];
  
  BD_LOOP({
    if(unif_rand() < 0.5 && nedgestype > 0) {
      // select an existing edge of type i at random, and propose toggling it off
      UnsrtELGetRand(Mtail, Mhead, sto->els[i]);
      
      // logratio is essentially copied from TNT, because the probability of 
      // choosing this particular mixing type cancels upon taking the ratio;
      // still need to count only edges and dyads of the appropriate mixing type, though
      MHp->logratio = log((nedgestype == 1 ? 1.0/(0.5*ndyadstype + 0.5) :
                           nedgestype / ((double) ndyadstype + nedgestype)));
    } else {
      // select a dyad of type i and propose toggling it        
      GetRandDyadFromLists(Mtail, // tail
                           Mhead, // head
                           sto->nodesbycode, // tails
                           sto->nodesbycode, // heads
                           sto->tailtypes + i, // tailattrs
                           sto->headtypes + i, // headattrs
                           sto->nodecountsbycode, // tailcounts
                           sto->nodecountsbycode, // headcounts
                           1, // length; only one allowed pairing since we've already sampled the strat type
                           ndyadstype, // dyadcount
                           TRUE, // diagonal; no higher level types here, so always TRUE
                           DIRECTED); // directed

      if(IS_OUTEDGE(Mtail[0],Mhead[0])) {
        // pick a new edge from the edgelist uniformly at random so we know its index
        // and hence don't have to look up the index of the edge tail -> head; this gives
        // the same probability of picking each existing edge as if we used the tail -> head
        // edge, but also allows us to keep the edgelists unsorted (at the cost of generating
        // an extra random index in this case)
        UnsrtELGetRand(Mtail, Mhead, sto->els[i]);

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
  UnsrtELToggleKnown(tail, head, sto->els[sto->currentmixingtype], edgeflag);
}

MH_F_FN(Mf_StratTNT) {
  // Free all the things
  GET_STORAGE(StratTNTStorage, sto);
  
  for(int i = 0; i < sto->nmixtypes; i++) {
    UnsrtELDestroy(sto->els[i]);
  }

  Free(sto->els);
  Free(sto->ndyadstype);
  Free(sto->nodesbycode);  
  
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


