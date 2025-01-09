/*  File inst/include/ergm_etamap.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#include "ergm_Rutil.h"

#ifndef _ERGM_ETAMAP_H_
#define _ERGM_ETAMAP_H_

void ergm_eta(double *theta, SEXP etamap, double *eta);
void ergm_etagrad(double *theta, SEXP etamap, double *eta);
void ergm_etagradmult(double *theta, double *v, unsigned int nv, SEXP etamap, double *ans);

#endif // _ERGM_ETAMAP_H_
