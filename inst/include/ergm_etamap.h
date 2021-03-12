#include "ergm_Rutil.h"

#ifndef _ERGM_ETAMAP_H_
#define _ERGM_ETAMAP_H_

void ergm_eta(double *theta, SEXP etamap, double *eta);
void ergm_etagrad(double *theta, SEXP etamap, double *eta);
void ergm_etagradmult(double *theta, double *v, unsigned int nv, SEXP etamap, double *ans);

#endif // _ERGM_ETAMAP_H_
