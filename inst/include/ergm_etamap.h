#include "ergm_Rutil.h"

void ergm_eta(double *theta, SEXP etamap, double *eta);
void ergm_etagrad(double *theta, SEXP etamap, double *eta);
void ergm_etagradmult(double *theta, double *v, unsigned int nv, SEXP etamap, double *ans);
