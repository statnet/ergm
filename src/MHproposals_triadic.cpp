/*  File src/MHproposals_triadic.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#define STRICT_MH_HEADERS

#include "MHproposals.h"
#include "ergm_dyad_hashmap.h"
#include "changestats_dgw_sp.h"

#include "cpp/ergm_network.h"
#include "cpp/ergm_changestat.h"
#include "cpp/ergm_proposal.h"

using ergm::ErgmCppNetwork;
using ergm::ErgmCppProposal;

/*********************
 void MH_SPDyad
*********************/
extern "C" MH_I_FN(Mi_SPDyad){
  Mi_TNT(MHp, nwp);
}

extern "C" MH_P_FN(Mp_SPDyad){
  ErgmCppNetwork nw(nwp);
  ErgmCppProposal<StoreDyadGenAndDegreeBoundAndModel> p(MHp);
  auto spcache = (StoreStrictDyadMapUInt *) p.aux_storage[0];

  // With probability 1-p.dinput[0], or if no dyad has any shared
  // partners, just fall back to TNT. This is OK to do because it is
  // impossible for the triadic proposal to produce a network with no
  // shared partners.
  if(unif_rand() > p.dinput[0]){
    Mp_TNT(MHp, nwp);
    return;
  }else if(kh_size(spcache) == 0){
    // It's triadic proposal's turn, but there isn't one to propose.
    p.tail[0] = MH_FAILED;
    p.head[0] = MH_CONSTRAINT;
    return;

    // FIXME: This code *should* adjust the acceptance probability if
    // the proposal falls back to TNT rather than failing, but it's
    // not working.

    /*   Rboolean edgestate = IS_OUTEDGE(*Mtail, p.head[0]); */
    /*   if(kh_size(spcache) == 0 && !edgestate){ */
    /*   // Find out if we are jumping from 0-SP to a 1-SP state: if we */
    /*   // currently have 0 SP and are adding an edge... */
    /*     Rboolean SP01 = FALSE; // Indicator of whether the toggle will create an SP of an appropriate type. */
    /*     switch(p.iinput[0]){ */
    /*     case L2UTP: */
    /*       SP01 = nw.degree[p.tail[0]] || nw.degree[p.head[0]]; break; */
    /*     case L2OTP: */
    /*     case L2ITP: */
    /*       SP01 = IN_DEG[*Mtail] || OUT_DEG[p.head[0]]; break; */
    /*     case L2OSP: SP01 = IN_DEG[p.head[0]] != 0; break; */
    /*     case L2ISP: SP01 = OUT_DEG[*Mtail] != 0; break; */
    /*     } */
    /*     if(!SP01) return; */

    /*     MHp->logratio += log1p(-p.dinput[0]) - log(DyadGenEdgecount(storage->gen)+1) */
    /*       - log(0.5) + log(storage->gen->ndyads); */
    /*     return; */
    /*   }else if(kh_size(spcache) == 1 && edgestate){ */
    /*   // Find out if we are jumping from 1-SP to a 0-SP state: if we */
    /*   // currently have 1 SP and are removing an edge... */
    /*     Rboolean SP10 = FALSE; // Indicator of whether the toggle will remove an SP of an appropriate type. */
    /*     switch(p.iinput[0]){ */
    /*     case L2UTP: */
    /*       SP10 = (nw.degree[p.tail[0]] == 1) || (nw.degree[p.head[0]] == 1); break; */
    /*     case L2OTP: */
    /*     case L2ITP: */
    /*       SP10 = (IN_DEG[*Mtail] == 1) || (OUT_DEG[p.head[0]] == 1); break; */
    /*     case L2OSP: SP10 = IN_DEG[p.head[0]] == 1; break; */
    /*     case L2ISP: SP10 = OUT_DEG[*Mtail] == 1; break; */
    /*     } */
    /*     if(!SP10) return; */

    /*     MHp->logratio += log(0.5) - log(storage->gen->ndyads) */
    /*       - log1p(-p.dinput[0]) + log(DyadGenEdgecount(storage->gen)); */
    /*     return; */
    /*   } */
    /* } */
  }

  L2Type type = (L2Type) p.iinput[0];

  BD_COND_LOOP(p.storage->bd, {
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
        if(type == L2ITP){
          p.tail[0] = dyad.head;
          p.head[0] = dyad.tail;
        }else{
          p.tail[0] = dyad.tail;
          p.head[0] = dyad.head;
        }
      }else{
        p.tail[0] = MIN(dyad.tail, dyad.head);
        p.head[0] = MAX(dyad.tail, dyad.head);
      }
    },
    DyadGenSearch(p.tail[0], p.head[0], p.storage->gen),
    4.0/MAX_TRIES);

  // If we keep trying to propose dyads that are fixed, fall back to
  // TNT.
  if(p.tail[0] == MH_FAILED && p.head[0] == MH_UNSUCCESSFUL){
    Mp_TNT(MHp, nwp);
    return;
  }

  CHECK_CHANGESTATS(p.storage->m);

  // q(y* | y) = 1/TD(y), where TD(y) is the number of transitive
  // dyads in y.
  Dyad oldtd = kh_size(spcache), newtd = oldtd;


  // The following is setting up to use macros developed for the *sp
  // terms.
  Rboolean edgeflag = nw(p.tail[0], p.head[0]);
  int echange = edgeflag ? -1 : +1;
  Vertex tail = p.tail[0], head = p.head[0];

#define sp_nonzero newtd += (L2 + echange != 0) - (L2 != 0);

  switch(type){
  case L2UTP: dspUTP_change(sp_nonzero, ); break;
  case L2OTP: dspOTP_change(sp_nonzero, ); break;
  case L2ITP: dspITP_change(sp_nonzero, ); break;
  case L2OSP: dspOSP_change(sp_nonzero, ); break;
  case L2ISP: dspISP_change(sp_nonzero, ); break;
  default: error("In ergm:Mp_SPDyad(), an unsupported type of triad: %d.", type);
  }

#undef sp_nonzero

  // q(y | y*) / q(y* | y) = 1/TD(y*) / (1/TD(y)) = TD(y) / TD(y*)
  p.logratio += log(oldtd) - log(newtd);
}

extern "C" MH_F_FN(Mf_SPDyad){
  Mf_TNT(MHp, nwp);
}
