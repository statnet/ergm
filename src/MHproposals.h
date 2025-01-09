/*  File src/MHproposals.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#ifndef MHPROPOSALS_H
#define MHPROPOSALS_H

#include "ergm_MHproposal.h"
#include "ergm_MHproposal_bd.h"
#include "ergm_dyadgen.h"

typedef struct{DyadGen *gen; DegreeBound *bd;} StoreDyadGenAndDegreeBound;

#define INIT_DYADGEN_AND_DEGREE_BOUND(el)                       \
  ALLOC_STORAGE(1, StoreDyadGenAndDegreeBound, storage);        \
  storage->gen = DyadGenInitializeR(MHp->R, nwp, el);           \
  storage->bd = DegreeBoundInitializeR(MHp->R, nwp);

#define DESTROY_DYADGEN_AND_DEGREE_BOUND                \
  GET_STORAGE(StoreDyadGenAndDegreeBound, storage);     \
  DyadGenDestroy(storage->gen);                         \
  DegreeBoundDestroy(storage->bd);

// Declared here so other routines can use it as a subroutine.
MH_I_FN(Mi_TNT);
MH_P_FN(Mp_TNT);
MH_F_FN(Mf_TNT);


#endif 
