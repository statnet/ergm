/*  File inst/include/inc/ergm_changestat.h.template.do_not_include_directly.h
 *  in package ergm, part of the Statnet suite of packages for network
 *  analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ETYPE(ModelTermstruct) {
  void (*c_func)(Vertex, Vertex, IFEWT(EWTTYPE,) struct ETYPE(ModelTermstruct)*, ETYPE(Network)*, EWTTYPE);
  void (*d_func)(Edge, Vertex*, Vertex*, IFEWT(EWTTYPE*,) struct ETYPE(ModelTermstruct)*, ETYPE(Network)*);
  void (*i_func)(struct ETYPE(ModelTermstruct)*, ETYPE(Network)*);
  void (*u_func)(Vertex, Vertex, IFEWT(EWTTYPE,) struct ETYPE(ModelTermstruct)*, ETYPE(Network)*, EWTTYPE);
  void (*f_func)(struct ETYPE(ModelTermstruct)*, ETYPE(Network)*);
  void (*s_func)(struct ETYPE(ModelTermstruct)*, ETYPE(Network)*);
  SEXP (*w_func)(struct ETYPE(ModelTermstruct)*, ETYPE(Network)*);  
  void (*x_func)(unsigned int, void *, struct ETYPE(ModelTermstruct)*, ETYPE(Network)*);
  void (*z_func)(struct ETYPE(ModelTermstruct)*, ETYPE(Network)*, Rboolean);
  double *attrib; /* Ptr to vector of covariates (if necessary; generally unused) */
  int *iattrib; /* Ptr to vector of integer covariates (if necessary; generally unused) */
  unsigned int nstats;   /* Number of change statistics to be returned */
  unsigned int statspos; /* Position of this term's stats in the workspace vector. */ 
  double *dstats; /* ptr to change statistics returned */
  unsigned int ninputparams; /* Number of double input parameters passed to function */
  double *inputparams; /* ptr to double input parameters passed */
  unsigned int niinputparams; /* Number of integer input parameters passed to function */
  int *iinputparams; /* ptr to integer input parameters passed */
  double *statcache; /* vector of the same length as dstats */
  double *emptynwstats; /* vector of the same length as dstats or NULL*/
  void *storage; /* optional space for persistent storage */
  void **aux_storage; /* optional space for persistent public (auxiliary) storage */
  unsigned int n_aux;
  unsigned int *aux_slots;
  SEXP R; /* R term object. */
  SEXP ext_state; /* A place from which to read extended state. */
} ETYPE(ModelTerm);

/****************************************************
 Macros to make life easier when writing C code for change statistics:  */

#ifdef __cplusplus
}
#endif

#include "ergm_changestat_common.do_not_include_directly.h"
