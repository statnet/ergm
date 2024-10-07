/*  File src/MHproposals_triadic.c in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2024 Statnet Commons
 */
#define STRICT_MH_HEADERS

#include "MHproposals.h"
#include "ergm_MHstorage.h"
#include "ergm_changestat.h"
#include "ergm_dyad_hashmap.h"
#include "changestats_dgw_sp.h"

/*********************
 void MH_SPDyad
*********************/
MH_I_FN(Mi_SPDyad){
  Mi_TNT(MHp, nwp);
}

MH_P_FN(Mp_SPDyad){
  MH_GET_STORAGE(StoreDyadGenAndDegreeBound, storage);
  MH_GET_AUX_STORAGE(StoreStrictDyadMapUInt, spcache);

  // With probability 1-MH_INPUTS[0], or if no dyad has any shared
  // partners, just fall back to TNT. This is OK to do because it is
  // impossible for the triadic proposal to produce a network with no
  // shared partners.
  if(unif_rand() > MH_INPUTS[0]){
    Mp_TNT(MHp, nwp);
    return;
  }else if(kh_size(spcache) == 0){
    // It's triadic proposal's turn, but there isn't one to propose.
    *Mtail = MH_FAILED;
    *Mhead = MH_CONSTRAINT;
    return;

    // FIXME: This code *should* adjust the acceptance probability if
    // the proposal falls back to TNT rather than failing, but it's
    // not working.

    /*   Rboolean edgestate = IS_OUTEDGE(*Mtail, *Mhead); */
    /*   if(kh_size(spcache) == 0 && !edgestate){ */
    /*   // Find out if we are jumping from 0-SP to a 1-SP state: if we */
    /*   // currently have 0 SP and are adding an edge... */
    /*     Rboolean SP01 = FALSE; // Indicator of whether the toggle will create an SP of an appropriate type. */
    /*     switch(MH_IINPUTS[0]){ */
    /*     case ESPUTP: */
    /*       SP01 = (OUT_DEG[*Mtail]+IN_DEG[*Mtail]) || (OUT_DEG[*Mhead]+IN_DEG[*Mhead]); break; */
    /*     case ESPOTP: */
    /*     case ESPITP: */
    /*       SP01 = IN_DEG[*Mtail] || OUT_DEG[*Mhead]; break; */
    /*     case ESPOSP: SP01 = IN_DEG[*Mhead] != 0; break; */
    /*     case ESPISP: SP01 = OUT_DEG[*Mtail] != 0; break; */
    /*     } */
    /*     if(!SP01) return; */

    /*     MHp->logratio += log1p(-MH_INPUTS[0]) - log(DyadGenEdgecount(storage->gen)+1) */
    /*       - log(0.5) + log(storage->gen->ndyads); */
    /*     return; */
    /*   }else if(kh_size(spcache) == 1 && edgestate){ */
    /*   // Find out if we are jumping from 1-SP to a 0-SP state: if we */
    /*   // currently have 1 SP and are removing an edge... */
    /*     Rboolean SP10 = FALSE; // Indicator of whether the toggle will remove an SP of an appropriate type. */
    /*     switch(MH_IINPUTS[0]){ */
    /*     case ESPUTP: */
    /*       SP10 = (OUT_DEG[*Mtail]+IN_DEG[*Mtail] == 1) || (OUT_DEG[*Mhead]+IN_DEG[*Mhead] == 1); break; */
    /*     case ESPOTP: */
    /*     case ESPITP: */
    /*       SP10 = (IN_DEG[*Mtail] == 1) || (OUT_DEG[*Mhead] == 1); break; */
    /*     case ESPOSP: SP10 = IN_DEG[*Mhead] == 1; break; */
    /*     case ESPISP: SP10 = OUT_DEG[*Mtail] == 1; break; */
    /*     } */
    /*     if(!SP10) return; */

    /*     MHp->logratio += log(0.5) - log(storage->gen->ndyads) */
    /*       - log1p(-MH_INPUTS[0]) + log(DyadGenEdgecount(storage->gen)); */
    /*     return; */
    /*   } */
    /* } */
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
      if(DIRECTED){
        if(MH_IINPUTS[0] == ESPITP){
          *Mtail = dyad.head;
          *Mhead = dyad.tail;
        }else{
          *Mtail = dyad.tail;
          *Mhead = dyad.head;
        }
      }else{
        *Mtail = MIN(dyad.tail, dyad.head);
        *Mhead = MAX(dyad.tail, dyad.head);
      }
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


  // The following is setting up to use macros developed for the *sp
  // terms.
  Rboolean edgeflag = IS_OUTEDGE(*Mtail, *Mhead);
  int echange = edgeflag ? -1 : +1;
  Vertex tail = *Mtail, head = *Mhead;

#define sp_nonzero newtd += (L + echange != 0) - (L != 0);

  switch(MH_IINPUTS[0]){
  case ESPUTP: dspUTP_change(L, sp_nonzero, ); break;
  case ESPOTP: dspOTP_change(L, sp_nonzero, ); break;
  case ESPITP: dspITP_change(L, sp_nonzero, ); break;
  case ESPOSP: dspOSP_change(L, sp_nonzero, ); break;
  case ESPISP: dspISP_change(L, sp_nonzero, ); break;
  }

#undef sp_nonzero

  // q(y | y*) / q(y* | y) = 1/TD(y*) / (1/TD(y)) = TD(y) / TD(y*)
  MHp->logratio += log(oldtd) - log(newtd);
}

MH_F_FN(Mf_SPDyad){
  Mf_TNT(MHp, nwp);
}
