/*  File src/MPLE.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2021 Statnet Commons
 */
#ifndef MPLE_H
#define MPLE_H

#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_model.h"
#include "ergm_rlebdm.h"
#include "ergm_state.h"

void MpleInit_hash_wl_RLE(ErgmState *s, int *responsevec, double *covmat, int *weightsvector,
			  RLEBDM1D *wl,
			  Edge maxNumDyads, Edge maxNumDyadTypes);
#endif
