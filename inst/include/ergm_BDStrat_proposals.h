#ifndef _ERGM_BDSTRAT_PROPOSALS_H_
#define _ERGM_BDSTRAT_PROPOSALS_H_

#include "ergm_MHproposal.h"
#include "ergm_edgetree_types.h"
#include "ergm_unsorted_edgelist.h"
#include "ergm_weighted_population.h"

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
  int *nodepos; // nodepos[i-1] is position of vertex i in nodesvec[vattr[i-1]]
  
  UnsrtEL *edgelist; // All edges in the network.
  
  int tailtype; // Attribute type of the last tail to be proposed.
  int tailmaxl; // Will the tail change the maximality status if the current proposal is accepted?
  
  int headtype; // Ditto for heads.
  int headmaxl;
  
  Dyad currentdyads; // Number of dyads that can be selected in the current network.
  Dyad proposeddyads; // As above, but if the proposal is accepted.
  
  int currentsubmaxledges;  // Number of edges in the current network both of whose endpoints are submaximal
  int proposedsubmaxledges; // Number of edges in the proposed network both of whose endpoints are submaximal
  
  int bound; // Single upper bound on degree.
  int nmixtypes; // Number of pairings of attributes.
  int *vattr; // Vertex attributes.
  // Parallel vectors of attribute combinations that are allowed.
  int *tailtypes;
  int *headtypes;
} BDTNTStorage;

typedef struct {
  UnsrtEL **els;
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
  
  int *BDtypesbyStrattype;
  int **BDtailsbyStrattype;
  int **BDheadsbyStrattype;
  
  int *strattailtypes;
  int *stratheadtypes;
  
  int nstratlevels;
  
  int *currentsubmaxledgestype;
  int **indmat;
  
  int nmixtypestoupdate;
  int *mixtypestoupdate;
} BDStratTNTStorage;



#endif 
