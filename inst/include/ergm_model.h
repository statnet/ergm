/*  File inst/include/ergm_model.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2019 Statnet Commons
 */
#ifndef _ERGM_MODEL_H_
#define _ERGM_MODEL_H_

#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "R_ext/Rdynload.h"
#include "ergm_MHproposal.h"

/* A Model object contains information about an entire ERGM, including the
   total numbers of terms, parameters, and statistics along with a pointer
   to an array of ModelTerm structures.  */
typedef struct Modelstruct {
  ModelTerm *termarray; /* array of size n_terms; see changestat.h
                           for ModelTerm definition */
  int n_terms;
  int n_stats;
  double *workspace; /* temporary workspace of size */
  double **dstatarray; /* array of size n_terms; the ith element in this
			  array is a pointer to an array of size
			  termarray[i].nstats                    */
  unsigned int n_aux;
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

#define UPDATE_STORAGE_COND(tail, head, nwp, m, MHp, cond){		\
    if(MHp && ((MHProposal*)MHp)->u_func) ((MHProposal*)MHp)->u_func(tail, head, MHp, nwp); \
    EXEC_THROUGH_TERMS(m, {						\
	IFDEBUG_BACKUP_DSTATS;						\
	if(mtp->u_func && (cond))					\
	  (*(mtp->u_func))(tail, head, mtp, nwp);  /* Call u_??? function */ \
	IFDEBUG_RESTORE_DSTATS;						\
      });								\
  }

#define UPDATE_STORAGE(tail, head, nwp, m, MHp){			\
    UPDATE_STORAGE_COND(tail, head, nwp, m, MHp, TRUE);			\
  }

Model* ModelInitialize (char *fnames, char *sonames, double **inputs,
			int n_terms);

void ModelDestroy(Network *nwp, Model *m);

/* A Model object contains information about an entire ERGM, including the
   total numbers of terms, parameters, and statistics along with a pointer
   to an array of ModelTerm structures.  */

int GetIndexForAttrValue(int value);

/* *** don't forget tail-> head, so this function accepts toggletail first, not togglehead  */

void ChangeStats(unsigned int ntoggles, Vertex *toggletail, Vertex *togglehead, Network *nwp, Model *m);

void InitStats(Network *nwp, Model *m);

void DestroyStats(Network *nwp, Model *m);

#endif

