#define STRICT_MH_HEADERS

#include "MHproposals.h"
#include "ergm_MHstorage.h"
#include "ergm_changestat.h"
#include "ergm_dyad_hashmap.h"

/*********************
 void MH_SPDyad
*********************/
MH_I_FN(Mi_SPDyad){
  Mi_TNT(MHp, nwp);
}

MH_P_FN(Mp_SPDyad){
  MH_GET_STORAGE(StoreDyadGenAndDegreeBound, storage);
  MH_GET_AUX_STORAGE(StoreDyadMapUInt, spcache);

  // With probability 1-MH_INPUTS[0], or if no dyad has any shared
  // partners, just fall back to TNT. This is OK to do because it is
  // impossible for the triadic proposal to produce a network with no
  // shared partners.
  if(kh_size(spcache) == 0 || unif_rand() > MH_INPUTS[0]){
    Mp_TNT(MHp, nwp);
    return;
  }

  BD_COND_LOOP(storage->bd, {
      // Select a random key from the shared partner hash table; only
      // those dyads that have at least one shared partner can be thus
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
    },
    DyadGenSearch(*Mtail, *Mhead, storage->gen),
    4.0/MAX_TRIES);

  // If we keep trying to propose dyads that are fixed, fall back to
  // TNT.
  if(*Mtail == MH_FAILED && *Mhead == MH_UNSUCCESSFUL){
    Mp_TNT(MHp, nwp);
    return;
  }

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

MH_F_FN(Mf_SPDyad){
  Mf_TNT(MHp, nwp);
}
