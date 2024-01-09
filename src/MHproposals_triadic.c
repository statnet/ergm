#define STRICT_MH_HEADERS

#include "Rmath.h"
#include "ergm_MHproposal.h"
#include "ergm_changestat.h"
#include "ergm_MHstorage.h"
#include "ergm_dyad_hashmap.h"


MH_I_FN(Mi_TNT_simple){
  MHp->ntoggles = 1;
}

MH_P_FN(Mp_TNT_simple){
  const double P=0.5, Q=1-P;
  double DP = P*N_DYADS, DO = DP/Q;

  Edge nedges = N_EDGES;
  double logratio=0;
  if (unif_rand() < P && nedges > 0) { /* Select a tie at random from the network of eligibles */
    GetRandEdge(Mtail, Mhead, nwp);
    logratio = TNT_LR_E(nedges, Q, DP, DO);
  }else{ /* Select a dyad at random from the list */
    GetRandDyad(Mtail, Mhead, nwp);

    if(IS_OUTEDGE(Mtail[0],Mhead[0])){
      logratio = TNT_LR_DE(nedges, Q, DP, DO);
    }else{
      logratio = TNT_LR_DN(nedges, Q, DP, DO);
    }
  }

  MHp->logratio += logratio;
}


/*********************
 void MH_SPDyad
*********************/
MH_I_FN(Mi_SPDyad){
  MHp->ntoggles = 1;
}

MH_P_FN(Mp_SPDyad){
  MH_GET_AUX_STORAGE(StoreDyadMapUInt, spcache);

  if(kh_size(spcache) == 0 || unif_rand() < 0.5){
    Mp_TNT_simple(MHp, nwp);
    return;
  }

  // Select a random key from the shared partner hash table; only
  // those dyads that have at least one shared partner can be
  // selected.
  khiter_t pos;
  do{
    pos = unif_rand() * kh_n_buckets(spcache);
  }while(!kh_exist(spcache, pos));
  TailHead dyad = kh_key(spcache, pos);

  // As of right now, the data structure does not guarantee correct
  // tail-head ordering for undirected networks.
  *Mtail = MIN(dyad.tail, dyad.head);
  *Mhead = MAX(dyad.tail, dyad.head);

  // q(y* | y) = 1/TD(y), where TD(y) is the number of transitive
  // dyads in y.
  Dyad oldtd = kh_size(spcache), newtd = oldtd;
  Rboolean edgeflag = IS_UNDIRECTED_EDGE(*Mtail, *Mhead);
  int echange = edgeflag ? -1 : +1;

  EXEC_THROUGH_EDGES(*Mhead,e,u, {
      if (u!=*Mtail){
        int L2tu = GETDMUI(*Mtail,u,spcache);
        newtd += (L2tu + echange != 0) - (L2tu != 0);
      }
    });

  EXEC_THROUGH_EDGES(*Mtail,e,u, {
      if (u!=*Mhead){
        int L2uh = GETDMUI(u,*Mhead,spcache);
        newtd += (L2uh + echange != 0) - (L2uh != 0);
      }
    });

  // q(y | y*) / q(y* | y) = 1/TD(y*) / (1/TD(y)) = TD(y) / TD(y*)
  MHp->logratio += log(oldtd) - log(newtd);
}

/* MH_F_FN(Mf_SPDyad){ */
/* } */
