/*  File inst/include/ergm_wtmodel.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2022 Statnet Commons
 */
#ifndef _ERGM_WTMODEL_H_
#define _ERGM_WTMODEL_H_

#include "ergm_wtedgetree.h"
#include "ergm_wtchangestat.h"
#include "R_ext/Rdynload.h"
#include "ergm_wtMHproposal.h"
#include "ergm_Rutil.h"

/* A WtModel object contains information about an entire ERGM, including the
   total numbers of terms, parameters, and statistics along with a pointer
   to an array of WtModelTerm structures.  */
typedef struct WtModelstruct {
  SEXP R; /* Pointer to the R ergm_model object. */
  SEXP ext_state; /* Pointer to the extended state for the whole model. */
  WtModelTerm *termarray; /* array of size n_terms; see changestat.h
                           for WtModelTerm definition */
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
} WtModel;

#define WtFOR_EACH_TERM(m) for(WtModelTerm *mtp = (m)->termarray; mtp < (m)->termarray + (m)->n_terms; mtp++)

#define WtEXEC_THROUGH_TERMS(m, subroutine)				\
  WtFOR_EACH_TERM(m){							\
    subroutine;								\
  }

#define WtFOR_EACH_TERM_INREVERSE(m) for(WtModelTerm *mtp = (m)->termarray + (m)->n_terms - 1; mtp >= (m)->termarray; mtp--)

#define WtEXEC_THROUGH_TERMS_INREVERSE(m, subroutine)			\
  WtFOR_EACH_TERM_INREVERSE(m){						\
    subroutine;								\
  }

#define WtEXEC_THROUGH_TERMS_INTO(m, output, subroutine)		\
  WtFOR_EACH_TERM(m){							\
    double *dstats = output + mtp->statspos;				\
    subroutine;								\
  }

#define WtSEND_X_SIGNAL(nwp, m, MHp, type, data)                        \
  if((MHp) && ((WtMHProposal*)(MHp))->x_func) ((WtMHProposal*)(MHp))->x_func((type), (data), (MHp), (nwp)); \
  WtEXEC_THROUGH_TERMS((m), {                                           \
      if(mtp->x_func)                                                   \
        (*(mtp->x_func))((type), (data), (mtp), (nwp));                 \
    });

#define WtSEND_X_SIGNAL_INREVERSE(nwp, m, MHp, type, data)              \
  WtEXEC_THROUGH_TERMS_INREVERSE((m), {                                 \
      if(mtp->x_func)                                                   \
        (*(mtp->x_func))((type), (data), (mtp), (nwp));                 \
    });                                                                 \
  if((MHp) && ((WtMHProposal*)(MHp))->x_func) ((WtMHProposal*)(MHp))->x_func((type), (data), (MHp), (nwp));

#define WtSEND_X_SIGNAL_INTO(nwp, m, MHp, output, type, data)           \
  if((MHp) && ((WtMHProposal*)(MHp))->x_func) ((WtMHProposal*)(MHp))->x_func((type), (data), (MHp), (nwp)); \
  WtEXEC_THROUGH_TERMS_INTO((m), output, {                              \
      if(mtp->x_func){                                                  \
        mtp->dstats = dstats;                                           \
        (*(mtp->x_func))((type), (data), (mtp), (nwp));                 \
      }                                                                 \
    });

 /* If DEBUG is set, back up mtp->dstats and set it to NULL in order
    to trigger a segfault if u_func tries to write to change
    statistics; then restore it. Otherwise, don't bother. */
#ifdef DEBUG
#define IFDEBUG_BACKUP_DSTATS double *dstats = mtp->dstats; mtp->dstats = NULL;
#define IFDEBUG_RESTORE_DSTATS mtp->dstats = dstats;
#else
#define IFDEBUG_BACKUP_DSTATS
#define IFDEBUG_RESTORE_DSTATS
#endif

WtModel* WtModelInitialize(SEXP mR, SEXP ext_stateR, WtNetwork *nwp, Rboolean noinit_s);

void WtModelDestroy(WtNetwork *nwp, WtModel *m);

/* A WtModel object contains information about an entire ERGM, including the
   total numbers of terms, parameters, and statistics along with a pointer
   to an array of WtModelTerm structures.  */

void WtChangeStats(unsigned int ntoggles, Vertex *tails, Vertex *heads, double *weights, WtNetwork *nwp, WtModel *m);
void WtChangeStats1(Vertex tail, Vertex head, double weight, WtNetwork *nwp, WtModel *m, double edgestate);
void WtZStats(WtNetwork *nwp, WtModel *m, Rboolean skip_s);
void WtEmptyNetworkStats(WtModel *m, Rboolean skip_s);
void WtSummStats(Edge n_edges, Vertex *tails, Vertex *heads, double *weights, WtNetwork *nwp, WtModel *m);
#endif

