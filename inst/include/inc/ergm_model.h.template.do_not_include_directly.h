/*  File inst/include/inc/ergm_model.h.template.do_not_include_directly.h in
 *  package ergm, part of the Statnet suite of packages for network analysis,
 *  https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#include "R_ext/Rdynload.h"
#include "../ergm_Rutil.h"

/* A ETYPE(Model) object contains information about an entire ERGM, including the
   total numbers of terms, parameters, and statistics along with a pointer
   to an array of ETYPE(ModelTerm) structures.  */
typedef struct ETYPE(Modelstruct) {
  SEXP R; /* Pointer to the R ergm_model object. */
  SEXP ext_state; /* Pointer to the extended state for the whole model. */
  ETYPE(ModelTerm) *termarray; /* array of size n_terms; see changestat.h
                           for ETYPE(ModelTerm) definition */
  int n_terms;
  int n_stats;
  unsigned int n_u; /* Number of terms with updaters. */
  double *workspace; /* temporary workspace of size n_stats */
  double *workspace_backup; /* since workspace is often replaced, we need to keep track of it for freeing */
  double **dstatarray; /* array of size n_terms; the ith element in this
			  array is a pointer to an array of size
			  termarray[i].nstats                    */
  unsigned int n_aux;
  Rboolean noinit_s;
} ETYPE(Model);

 /* If NDEBUG is unset, back up mtp->dstats and set it to NULL in
    order to trigger a segfault if u_func tries to write to change
    statistics; then restore it. Otherwise, don't bother. */
#ifndef NDEBUG
#define IFDEBUG_BACKUP_DSTATS double *dstats = mtp->dstats; mtp->dstats = NULL;
#define IFDEBUG_RESTORE_DSTATS mtp->dstats = dstats;
#else
#define IFDEBUG_BACKUP_DSTATS
#define IFDEBUG_RESTORE_DSTATS
#endif

ETYPE(Model)* ETYPE(ModelInitialize)(SEXP mR, SEXP ext_stateR, ETYPE(Network) *nwp, Rboolean noinit_s);

void ETYPE(ModelDestroy)(ETYPE(Network) *nwp, ETYPE(Model) *m);

/* A ETYPE(Model) object contains information about an entire ERGM, including the
   total numbers of terms, parameters, and statistics along with a pointer
   to an array of ETYPE(ModelTerm) structures.  */

void ETYPE(ChangeStatsDo)(unsigned int ntoggles, Vertex *tails, Vertex *heads, IFEWT(EWTTYPE *weights,) ETYPE(Network) *nwp, ETYPE(Model) *m);
void ETYPE(ChangeStatsUndo)(unsigned int ntoggles, Vertex *tails, Vertex *heads, IFEWT(EWTTYPE *weights,) ETYPE(Network) *nwp, ETYPE(Model) *m);
void ETYPE(ChangeStats)(unsigned int ntoggles, Vertex *tails, Vertex *heads, IFEWT(EWTTYPE *weights,) ETYPE(Network) *nwp, ETYPE(Model) *m);
void ETYPE(ChangeStats1)(Vertex tail, Vertex head, IFEWT(EWTTYPE weight,) ETYPE(Network) *nwp, ETYPE(Model) *m, EWTTYPE edgestate);
void ETYPE(ZStats)(ETYPE(Network) *nwp, ETYPE(Model) *m, Rboolean skip_s);
void ETYPE(EmptyNetworkStats)(ETYPE(Model) *m, Rboolean skip_s);
void ETYPE(SummStats)(Edge n_edges, Vertex *tails, Vertex *heads, IFEWT(EWTTYPE *weights,) ETYPE(Network) *nwp, ETYPE(Model) *m);

