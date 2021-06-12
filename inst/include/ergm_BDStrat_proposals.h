/*  File inst/include/ergm_BDStrat_proposals.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2021 Statnet Commons
 */
#ifndef _ERGM_BDSTRAT_PROPOSALS_H_
#define _ERGM_BDSTRAT_PROPOSALS_H_

#include "ergm_MHproposal.h"
#include "ergm_edgetree_types.h"
#include "ergm_unsorted_edgelist.h"
#include "ergm_weighted_population.h"
#include "ergm_hash_edgelist.h"
#include "ergm_BDStratBlocks.h"

#define OUTVAL_NET(e,n) ((n)->outedges[(e)].value)
#define INVAL_NET(e,n) ((n)->inedges[(e)].value)
#define MIN_OUTEDGE_NET(a,n) (EdgetreeMinimum((n)->outedges, (a)))
#define MIN_INEDGE_NET(a,n) (EdgetreeMinimum((n)->inedges, (a)))
#define NEXT_OUTEDGE_NET(e,n) (EdgetreeSuccessor((n)->outedges,(e)))
#define NEXT_INEDGE_NET(e,n) (EdgetreeSuccessor((n)->inedges,(e)))
#define STEP_THROUGH_OUTEDGES_NET(a,e,v,n) for((e)=MIN_OUTEDGE_NET((a),(n));((v)=OUTVAL_NET((e),(n)))!=0;(e)=NEXT_OUTEDGE_NET((e),(n)))
#define STEP_THROUGH_INEDGES_NET(a,e,v,n) for((e)=MIN_INEDGE_NET((a),(n));((v)=INVAL_NET((e),(n)))!=0;(e)=NEXT_INEDGE_NET((e),(n)))

#define EXEC_THROUGH_EDGES_EA_NET_DECL(node, ego, alter, edge, net, subroutine) { \
  Vertex (ego) = (node); \
  Vertex (alter); \
  Edge (edge); \
  for((edge) = MIN_OUTEDGE_NET((ego),(net)); ((alter) = OUTVAL_NET((edge),(net))) != 0; (edge) = NEXT_OUTEDGE_NET((edge),(net))) { \
    subroutine \
  } \
  for((edge) = MIN_INEDGE_NET((ego),(net)); ((alter) = INVAL_NET((edge),(net))) != 0; (edge) = NEXT_INEDGE_NET((edge),(net))) { \
    subroutine \
  } \
}

#define EXEC_THROUGH_EDGES_EATH_NET_DECL(node, ego, alter, tail, head, edge, net, subroutine) { \
  Vertex (ego) = (node); \
  Vertex (alter), (tail), (head); \
  Edge (edge); \
  (tail) = (ego); \
  for((edge) = MIN_OUTEDGE_NET((ego),(net)); ((head) = (alter) = OUTVAL_NET((edge),(net))) != 0; (edge) = NEXT_OUTEDGE_NET((edge),(net))) { \
    subroutine \
  } \
  (head) = (ego); \
  for((edge) = MIN_INEDGE_NET((ego),(net)); ((tail) = (alter) = INVAL_NET((edge),(net))) != 0; (edge) = NEXT_INEDGE_NET((edge),(net))) { \
    subroutine \
  } \
}

MH_I_FN(Mi_BDStratTNT);

MH_F_FN(Mf_BDStratTNT);

typedef struct {
  BDStratBlocks *blocks;

  HashEL **hash;
  
  int tailmaxl;  
  int headmaxl;
  
  int stratmixingtype;
  
  double currentcumprob;
  double proposedcumprob;
  
  double *originalprobvec;
  
  WtPop *wtp;
  
  int **maxout;
  int **maxin;
  
  int **indegree;
  int **outdegree;
  int nbdlevels;
  
  int nmixtypes;
  
  int *strat_vattr;
  int *blocks_vattr;
  int *bd_vattr;  
  
  int nstratlevels;
  
  int **indmat;
    
  int nmixtypestoupdate;
  int *mixtypestoupdate;

  int CD;
} BDStratTNTStorage;

// determines which strat mixing types (if any) will have a change in toggleability status if we make the proposed toggle
static inline void ComputeChangesToToggleability(Vertex *tail, Vertex *head, BDStratTNTStorage *sto) {
  // here we compute the proposedcumprob, checking only those
  // mixing types that can be influenced by toggles made on 
  // the current mixing type
  sto->proposedcumprob = sto->currentcumprob;
  sto->nmixtypestoupdate = 0; // reset counter
  // avoid these somewhat expensive checks in the typical case
  // where you have enough submaximal nodes that you cannot
  // be exhausting any mixing types of toggleable dyads
  if(sto->blocks->last_tails->length <= 2 || sto->blocks->last_heads->length <= 2) {
    // temporarily set tail and head toggleability to what it would be in the proposed network
    BDStratBlocksToggleIf(*tail, *head, sto->blocks, sto->tailmaxl, sto->headmaxl);
    
    // how many strat types do we need to check?
    int ntocheck = ((!sto->blocks->directed) && (sto->strat_vattr[*tail] == sto->strat_vattr[*head])) ? sto->nstratlevels : 2*sto->nstratlevels;

    for(int i = 0; i < ntocheck; i++) {
      // find the index of the i'th strat type we need to check, by looking it up in the indmat
      int infl_i = sto->indmat[i < sto->nstratlevels ? sto->strat_vattr[*tail] : i - sto->nstratlevels][i < sto->nstratlevels ? i : sto->strat_vattr[*head]];

      // if this strat type is not included in the proposal, or is the same as the strat type of the proposed toggle,
      // then it cannot change toggleability status, so skip it
      if(infl_i < 0 || infl_i == sto->stratmixingtype) {
        continue;
      }
      
      // can we toggle this mixing type in the current network?
      int toggle_curr = WtPopGetWt(infl_i, sto->wtp) > 0;
      
      // will we be able to toggle this mixing type in the proposed network? 
      int toggle_prop = sto->hash[infl_i]->list->nedges > 0 || BDStratBlocksDyadCountPositive(sto->blocks, infl_i);
      
      // will there be a change in toggleability status?
      int change = toggle_curr - toggle_prop;

      // if so, take this into account      
      if(change) {
        sto->proposedcumprob -= change*sto->originalprobvec[infl_i];
        sto->mixtypestoupdate[sto->nmixtypestoupdate] = infl_i;
        sto->nmixtypestoupdate++;        
      }
    }
    
    // restore tail and head toggleability to their current status
    BDStratBlocksToggleIf(*tail, *head, sto->blocks, sto->tailmaxl, sto->headmaxl);
  }
}

#endif 
