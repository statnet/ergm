/*  File src/wtmodel.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2013 Statnet Commons
 */
#ifndef WTMODEL_H
#define WTMODEL_H

#include "wtedgetree.h"
#include "wtchangestat.h"
#include "R_ext/Rdynload.h"

/* A WtModel object contains information about an entire ERGM, including the
   total numbers of terms, parameters, and statistics along with a pointer
   to an array of WtModelTerm structures.  */
typedef struct WtModelstruct {
  WtModelTerm *termarray; /* array of size n_terms; see changestat.h
                           for WtModelTerm definition */
  int n_terms;
  int n_stats;
  double *workspace; /* temporary workspace of size */
  double **dstatarray; /* array of size n_terms; the ith element in this
			  array is a pointer to an array of size
			  termarray[i].nstats                    */
  unsigned int n_aux;
} WtModel;

#define WtFOR_EACH_TERM for(WtModelTerm *mtp = m->termarray; mtp < m->termarray + m->n_terms; mtp++)

#define WtEXEC_THROUGH_TERMS(subroutine){				\
    WtFOR_EACH_TERM{							\
      subroutine;							\
    }									\
  }

#define WtFOR_EACH_TERM_INREVERSE for(WtModelTerm *mtp = m->termarray + m->n_terms - 1; mtp >= m->termarray; mtp--)

#define WtEXEC_THROUGH_TERMS_INREVERSE(subroutine){			\
    WtFOR_EACH_TERM_INREVERSE{						\
      subroutine;							\
    }									\
  }


#define WtEXEC_THROUGH_TERMS_INTO(output, subroutine){			\
    double *dstats = output;						\
    WtFOR_EACH_TERM{							\
      subroutine;							\
      dstats += mtp->nstats;						\
    }									\
  }

 /* If DEBUG is set, back up mtp->dstats and set it to NULL in order
    to trigger a segfault if u_func tries to write to change
    statistics; then restore it. Otherwise, don't bother. */
#ifdef DEBUG

#define WtUPDATE_STORAGE_COND(tail, head, weight, m, nwp, cond){	\
    WtEXEC_THROUGH_TERMS({						\
	double *dstats = mtp->dstats; /* Back up mtp->dstats. */	\
	mtp->dstats = NULL; /* Trigger segfault if u_func tries to write to change statistics. */ \
	if(mtp->u_func && (cond))					\
	  (*(mtp->u_func))(tail, head, weight, mtp, nwp);  /* Call u_??? function */ \
	mtp->dstats = dstats; /* Restore mtp->dstats. */		\
      });								\
  }

#else


#define WtUPDATE_STORAGE_COND(tail, head, weight, m, nwp, cond){	\
    WtEXEC_THROUGH_TERMS({						\
	if(mtp->u_func && (cond))					\
	  (*(mtp->u_func))(tail, head, weight, mtp, nwp);  /* Call u_??? function */ \
      });								\
  }

#endif

#define WtUPDATE_STORAGE(tail, head, weight, m, nwp){			\
    WtUPDATE_STORAGE_COND(tail, head, weight, m, nwp, TRUE);		\
  }

WtModel* WtModelInitialize (char *fnames, char *sonames, double **inputs,
			int n_terms);

void WtModelDestroy(WtModel *m, WtNetwork *nwp);

/* A WtModel object contains information about an entire ERGM, including the
   total numbers of terms, parameters, and statistics along with a pointer
   to an array of WtModelTerm structures.  */

void WtChangeStats(unsigned int ntoggles, Vertex *toggletail, Vertex *togglehead, double *toggleweight, WtNetwork *nwp, WtModel *m);

void WtInitStats(WtNetwork *nwp, WtModel *m);

void WtDestroyStats(WtNetwork *nwp, WtModel *m);



#endif

