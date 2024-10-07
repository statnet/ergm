/*  File inst/include/ergm_model.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2024 Statnet Commons
 */
#ifndef _ERGM_MODEL_H_
#define _ERGM_MODEL_H_

#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "R_ext/Rdynload.h"
#include "ergm_MHproposal.h"
#include "ergm_Rutil.h"

/* A Model object contains information about an entire ERGM, including the
   total numbers of terms, parameters, and statistics along with a pointer
   to an array of ModelTerm structures.  */
typedef struct Modelstruct {
  SEXP R; /* Pointer to the R ergm_model object. */
  SEXP ext_state; /* Pointer to the extended state for the whole model. */
  ModelTerm *termarray; /* array of size n_terms; see changestat.h
                           for ModelTerm definition */
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
} Model;

#define FOR_EACH_TERM(m) for(ModelTerm *mtp = (m)->termarray; mtp < (m)->termarray + (m)->n_terms; mtp++)


#define EXEC_THROUGH_TERMS(m, subroutine)				\
  FOR_EACH_TERM(m){							\
    subroutine;								\
  }

#define FOR_EACH_TERM_INREVERSE(m) for(ModelTerm *mtp = (m)->termarray + (m)->n_terms - 1; mtp >= (m)->termarray; mtp--)

#define EXEC_THROUGH_TERMS_INREVERSE(m, subroutine)			\
  FOR_EACH_TERM_INREVERSE(m){						\
    subroutine;								\
  }

#define EXEC_THROUGH_TERMS_INTO(m, output, subroutine)			\
  FOR_EACH_TERM(m){							\
    double *dstats = output + mtp->statspos;				\
    subroutine;								\
  }

#define SEND_X_SIGNAL(nwp, m, MHp, type, data)                          \
  if((MHp) && ((MHProposal*)(MHp))->x_func) ((MHProposal*)(MHp))->x_func((type), (data), (MHp), (nwp)); \
  EXEC_THROUGH_TERMS((m), {                                             \
      if(mtp->x_func)                                                   \
        (*(mtp->x_func))((type), (data), (mtp), (nwp));                 \
    });

#define SEND_X_SIGNAL_INREVERSE(nwp, m, MHp, type, data)                 \
  EXEC_THROUGH_TERMS_INREVERSE((m), {                                   \
      if(mtp->x_func)                                                   \
        (*(mtp->x_func))((type), (data), (mtp), (nwp));                 \
    });                                                                 \
  if((MHp) && ((MHProposal*)(MHp))->x_func) ((MHProposal*)(MHp))->x_func((type), (data), (MHp), (nwp));

#define SEND_X_SIGNAL_INTO(nwp, m, MHp, output, type, data)              \
  if((MHp) && ((MHProposal*)(MHp))->x_func) ((MHProposal*)(MHp))->x_func((type), (data), (MHp), (nwp)); \
  EXEC_THROUGH_TERMS_INTO((m), output, {                                \
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

Model* ModelInitialize(SEXP mR, SEXP ext_stateR, Network *nwp, Rboolean noinit_s);

void ModelDestroy(Network *nwp, Model *m);

/* A Model object contains information about an entire ERGM, including the
   total numbers of terms, parameters, and statistics along with a pointer
   to an array of ModelTerm structures.  */

int GetIndexForAttrValue(int value);

/* *** don't forget tail-> head, so this function accepts toggletail first, not togglehead  */

void ChangeStats(unsigned int ntoggles, Vertex *tails, Vertex *heads, Network *nwp, Model *m);
void ChangeStats1(Vertex tail, Vertex head, Network *nwp, Model *m, Rboolean edgestate);
void ZStats(Network *nwp, Model *m, Rboolean skip_s);
void EmptyNetworkStats(Model *m, Rboolean skip_s);
void SummStats(Edge n_edges, Vertex *tails, Vertex *heads, Network *nwp, Model *m);
#endif

