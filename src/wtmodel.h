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

#define UPDATE_STORAGE(tail, head, weight, m, nwp){			\
  WtModelTerm *mtp = m->termarray;					\
  for (unsigned int i=0; i < m->n_terms; i++, mtp++){			\
  mtp->dstats = NULL; /* Trigger segfault if u_func tries to write to change statistics. */ \
  if(mtp->u_func) /* Skip if either no update or it's a d_func, so it doesn't require storage updates for provisional updates. */ \
    (*(mtp->u_func))(tail, head, weight, mtp, nwp);  /* Call u_??? function */ \
  }									\
  }

#define UPDATE_C_STORAGE(tail, head, weight, m, nwp){			\
  WtModelTerm *mtp = m->termarray;					\
  for (unsigned int i=0; i < m->n_terms; i++, mtp++){			\
  mtp->dstats = NULL; /* Trigger segfault if u_func tries to write to change statistics. */ \
  if(mtp->u_func && mtp->d_func==NULL) /* Skip if either no update or it's a d_func, so it doesn't require storage updates for provisional updates. */ \
    (*(mtp->u_func))(tail, head, weight, mtp, nwp);  /* Call u_??? function */ \
  }									\
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

