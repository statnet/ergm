#ifndef _ERGM_BDSTRAT_PROPOSALS_H_
#define _ERGM_BDSTRAT_PROPOSALS_H_

#include "ergm_MHproposal.h"
#include "ergm_edgetree_types.h"
#include "ergm_unsorted_edgelist.h"
#include "ergm_weighted_population.h"
#include "ergm_hash_edgelist.h"

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
  UnsrtEL **els;
  int currentmixingtype;
  Vertex **nodesbycode;
  
  int nmixtypes;

  WtPop *wtp;
  
  int *tailtypes;
  int *headtypes;

  Dyad *ndyadstype;
  int *nodecountsbycode;
} StratTNTStorage;

typedef struct {
  int *attrcounts; // Count of the number of nodes with each attribute type i that are "submaximal degree" (attrcounts[i] lengths of nodesvec[i]).
  Vertex **nodesvec; // List of lists of submaximal nodes of attribute i.
  int *nodepos; // nodepos[i] is position of vertex i in nodesvec[vattr[i]]
  
  HashEL *hash; // All edges in the network.
  
  int tailtype; // Attribute type of the last tail to be proposed.
  int tailmaxl; // Will the tail change the maximality status if the current proposal is accepted?
  
  int headtype; // Ditto for heads.
  int headmaxl;
  
  Dyad currentdyads; // Number of dyads that can be selected in the current network.
  Dyad proposeddyads; // As above, but if the proposal is accepted.
    
  int bound; // Single upper bound on degree.
  int nmixtypes; // Number of pairings of attributes.
  int *vattr; // Vertex attributes.
  int nlevels; // number of attribute levels
  
  // Parallel vectors of attribute combinations that are allowed.
  int *tailtypes;
  int *headtypes;
} BDTNTStorage;

typedef struct {
  HashEL **hash;
  Vertex ***nodesvec;
  int **attrcounts;
  
  int strattailtype;
  int bdtailtype;
  int tailmaxl;
  
  int stratheadtype;
  int bdheadtype;  
  int headmaxl;

  int *nodepos;
  
  int stratmixingtype;
  
  double currentcumprob;
  double proposedcumprob;
  
  double *originalprobvec;
  
  WtPop *wtp;
  
  int bound;
  int nmixtypes;
  
  int *strat_vattr;
  int *bd_vattr;
  
  int *bd_mixtypes;
  int *bd_tails;
  int *bd_heads;
  
  int *strattailtypes;
  int *stratheadtypes;
  
  int nstratlevels;
  int nbdlevels;
  
  int **indmat;
    
  int nmixtypestoupdate;
  int *mixtypestoupdate;
} BDStratTNTStorage;



#endif 
