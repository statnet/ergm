/*  File src/changestats.h in package ergm, part of the Statnet suite of
 *  packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */
#ifndef CHANGESTATS_H
#define CHANGESTATS_H

#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_Rutil.h"

static inline void cutoff_error(ModelTerm *mtp){
  error("%s", CHAR(STRING_ELT(getListElement(mtp->R, "cutoff.message"), 0)));
}

#endif
