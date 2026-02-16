/*  File inst/include/ergm_BDStrat_proposals.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */
#ifndef _ERGM_BDSTRAT_PROPOSALS_H_
#define _ERGM_BDSTRAT_PROPOSALS_H_

#include "ergm_MHproposal.h"
#include "ergm_edgetree_types.h"
#include "ergm_unsorted_edgelist.h"
#include "ergm_weighted_population.h"
#include "ergm_hash_edgelist.h"
#include "ergm_BDStratBlocks.h"

// initialization and finalization functions
// (also used by the temporal extension of this proposal)
MH_I_FN(Mi_BDStratTNT);
MH_F_FN(Mf_BDStratTNT);

// overall data structure for BDStratTNT proposal
typedef struct {
  // submaximal nodes
  BDNodeLists *lists;

  // submaximal dyads
  BDStratBlocks *blocks;

  // edges
  HashEL **hash;

  // degree bounds
  int **maxout;
  int **maxin;

  // current degrees
  int **indegree;
  int **outdegree;

  // current strat mixing type
  int stratmixingtype;

  // will tail/head change maximality on toggle?
  int tailmaxl;
  int headmaxl;

  // count and indices of strat mixing types whose toggleability
  // will change if we accept the current proposal
  int strat_nmixtypestoupdate;
  int *strat_mixtypestoupdate;

  // total weight of toggleable strat mixing types
  // in the current and proposed networks
  double current_total_weight;
  double proposed_total_weight;

  // original weights for strat mixing types
  double *original_weights;

  // weighted sampling data structure for strat mixing types
  WtPop *wtp;

  // vertex attribute values
  int *strat_vattr;
  int *blocks_vattr;
  int *bd_vattr;

  // basic attribute information (number of levels, number of included mixing types)
  int bd_nlevels;
  int strat_nlevels;
  int strat_nmixtypes;

  // matrix storing indices of included strat mixing types
  int **indmat;

  // flag: are we running the contrastive divergence algorithm?
  int CD;
} BDStratTNTStorage;

// helper function to determine which strat mixing types (if any) will
// have a change in toggleability status if we accept the proposed toggle
static inline void ComputeChangesToToggleability(Vertex *tail, Vertex *head, BDStratTNTStorage *sto) {
  // here we compute the proposed_total_weight, checking only those
  // mixing types that can be influenced by toggles made on 
  // the current mixing type
  sto->proposed_total_weight = sto->current_total_weight;
  sto->strat_nmixtypestoupdate = 0; // reset counter
  // avoid these somewhat expensive checks in the typical case
  // where you have enough submaximal nodes that you cannot
  // be exhausting any mixing types of toggleable dyads
  int ntails = BDNodeListsTailCount(*tail, *head, sto->lists);
  int nheads = BDNodeListsHeadCount(*tail, *head, sto->lists);
  if(ntails <= 2 || nheads <= 2) {
    // temporarily set tail and head toggleability to what it
    // would be in the proposed network
    BDNodeListsToggleIf(*tail, *head, sto->lists, sto->tailmaxl, sto->headmaxl);

    // how many strat types do we need to check?
    int tailattr = sto->strat_vattr[*tail];
    int headattr = sto->strat_vattr[*head];
    int ntocheck = ((!sto->lists->directed) && (tailattr == headattr)) 
                   ? sto->strat_nlevels : 2*sto->strat_nlevels;

    for(int i = 0; i < ntocheck; i++) {
      // find the index of the i'th strat type we need to check,
      // by looking it up in the indmat
      int row = i < sto->strat_nlevels ? tailattr : i - sto->strat_nlevels;
      int col = i < sto->strat_nlevels ? i : headattr;
      int infl_i = sto->indmat[row][col];

      // if this strat type is not included in the proposal, or is the same as
      // the strat type of the proposed toggle, then it cannot change
      // toggleability status, so skip it
      if(infl_i < 0 || infl_i == sto->stratmixingtype) {
        continue;
      }

      // can we toggle this mixing type in the current network?
      int toggle_curr = WtPopGetWt(infl_i, sto->wtp) > 0;

      // will we be able to toggle this mixing type in the proposed network? 
      int toggle_prop = HashELSize(sto->hash[infl_i]) > 0
                        || BDStratBlocksDyadCountPositive(sto->blocks, infl_i);

      // will there be a change in toggleability status?
      int change = toggle_prop - toggle_curr;

      // if so, take this into account
      if(change) {
        sto->proposed_total_weight += change*sto->original_weights[infl_i];
        sto->strat_mixtypestoupdate[sto->strat_nmixtypestoupdate] = infl_i;
        sto->strat_nmixtypestoupdate++;
      }
    }

    // restore tail and head toggleability to their current status
    BDNodeListsToggleIf(*tail, *head, sto->lists, sto->tailmaxl, sto->headmaxl);
  }
}

#endif
