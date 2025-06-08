/*  File src/ergm_Rutil.c in package ergm, part of the Statnet suite of
 *  packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#include <R.h>
#include "ergm_Rutil.h"
#include "ergm_constants.h"

int GetBuiltABIVersion_ergm(void){
  return ERGM_ABI_VERSION;
}

SEXP GetBuiltABIVersion_wrapper(SEXP client, SEXP lib){
  const char *sn = FIRSTCHAR(client), *ln = FIRSTCHAR(lib);
  int lfn = strlen("GetBuiltABIVersion_") + strlen(ln);
  char *fn = R_Calloc(lfn + 1, char);

  snprintf(fn, lfn + 1, "GetBuiltABIVersion_%s", ln);
  for(unsigned int i = 0; i < lfn; i++) if(fn[i] == '.') fn[i] = '_';

  int (*f)(void) = (int (*)(void)) R_FindSymbol(fn, sn, NULL);

  R_Free(fn);

  return f ?  ScalarInteger(f()) : R_NilValue;
}
