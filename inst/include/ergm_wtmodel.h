/*  File inst/include/ergm_wtmodel.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2020 Statnet Commons
 */
#ifndef WTMODEL_H
#define WTMODEL_H

#include "ergm_wtedgetree.h"
#include "ergm_wtchangestat.h"
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
} WtModel;

WtModel* WtModelInitialize (char *fnames, char *sonames, double **inputs,
			int n_terms);

void WtModelDestroy(WtModel *m);

/* A WtModel object contains information about an entire ERGM, including the
   total numbers of terms, parameters, and statistics along with a pointer
   to an array of WtModelTerm structures.  */

void WtChangeStats(unsigned int ntoggles, Vertex *toggletail, Vertex *togglehead, double *toggleweight, WtNetwork *nwp, WtModel *m);

#endif

