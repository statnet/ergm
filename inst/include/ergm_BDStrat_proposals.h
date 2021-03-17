#ifndef _ERGM_BDSTRAT_PROPOSALS_H_
#define _ERGM_BDSTRAT_PROPOSALS_H_

#include "ergm_MHproposal.h"
#include "ergm_edgetree_types.h"
#include "ergm_unsorted_edgelist.h"
#include "ergm_weighted_population.h"
#include "ergm_hash_edgelist.h"
#include "ergm_nodelist.h"

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

MH_I_FN(Mi_StratTNT);
MH_I_FN(Mi_BDTNT);
MH_I_FN(Mi_BDStratTNT);

MH_F_FN(Mf_StratTNT);
MH_F_FN(Mf_BDTNT);
MH_F_FN(Mf_BDStratTNT);

typedef struct {
  NodeList *nodelist;
  
  UnsrtEL **els;

  WtPop *wtp;
  
  Dyad *ndyadstype;

  int **indmat;
  
  int *vattr;

  int nmixtypes;
  int currentmixingtype;
  
  int CD;
} StratTNTStorage;

typedef struct {
  NodeList *nodelist;
  
  HashEL *hash; // All edges in the network.
  
  int tailmaxl; // Will the tail change the maximality status if the current proposal is accepted?  
  int headmaxl; // ditto head
  
  Dyad currentdyads; // Number of dyads that can be selected in the current network.
  Dyad proposeddyads; // As above, but if the proposal is accepted.
    
  int bound; // Single upper bound on degree.

  int CD;
} BDTNTStorage;

typedef struct {
  NodeList *nodelist;

  HashEL **hash;
  
  int tailmaxl;  
  int headmaxl;
  
  int stratmixingtype;
  
  double currentcumprob;
  double proposedcumprob;
  
  double *originalprobvec;
  
  WtPop *wtp;
  
  int bound;
  int nmixtypes;
  
  int *strat_vattr;
  int *bd_vattr;
  
  int nstratlevels;
  
  int **indmat;
    
  int nmixtypestoupdate;
  int *mixtypestoupdate;

  int CD;
} BDStratTNTStorage;

// determines which strat mixing types (if any) will have a change in toggleability status if we make the proposed toggle
static inline void ComputeChangesToToggleability(Vertex *tail, Vertex *head, int edgeflag, BDStratTNTStorage *sto) {
  // here we compute the proposedcumprob, checking only those
  // mixing types that can be influenced by toggles made on 
  // the current mixing type
  sto->proposedcumprob = sto->currentcumprob;
  sto->nmixtypestoupdate = 0; // reset counter
  // avoid these somewhat expensive checks in the typical case
  // where you have enough submaximal nodes that you cannot
  // be exhausting any mixing types of toggleable dyads
  if(sto->nodelist->attrcounts[sto->strat_vattr[*tail]][sto->bd_vattr[*tail]] <= 2 || sto->nodelist->attrcounts[sto->strat_vattr[*head]][sto->bd_vattr[*head]] <= 2) {
    // temporarily set tail and head toggleability to what it would be in the proposed network
    NodeListToggleKnownIf(*tail, *head, sto->nodelist, !edgeflag, sto->tailmaxl, sto->headmaxl);
    
    // how many strat types do we need to check?
    int ntocheck = (sto->strat_vattr[*tail] == sto->strat_vattr[*head]) ? sto->nstratlevels : 2*sto->nstratlevels;

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
      int toggle_prop = sto->hash[infl_i]->list->nedges > 0 || NodeListDyadCountPositive(sto->nodelist, infl_i);
      
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
    NodeListToggleKnownIf(*tail, *head, sto->nodelist, edgeflag, sto->tailmaxl, sto->headmaxl);
  }
}

#endif 
