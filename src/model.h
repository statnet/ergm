/*
 *  File ergm/src/model.h
 *  Part of the statnet package, http://statnetproject.org
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) in
 *    http://statnetproject.org/attribution
 *
 * Copyright 2003 Mark S. Handcock, University of Washington
 *                David R. Hunter, Penn State University
 *                Carter T. Butts, University of California - Irvine
 *                Steven M. Goodreau, University of Washington
 *                Martina Morris, University of Washington
 * Copyright 2007 The statnet Development Team
 */
#ifndef MODEL_H
#define MODEL_H

#include "edgetree.h"
#include "changestats.h"
#include "R_ext/Rdynload.h"

/* A Model object contains information about an entire ERGM, including the
   total numbers of terms, parameters, and statistics along with a pointer
   to an array of ModelTerm structures.  */
typedef struct Modelstruct {
  ModelTerm *termarray; /* array of size n_terms; see changestats.h
                           for ModelTerm definition */
  int n_terms;
  int n_stats;
  double *workspace; /* temporary workspace of size */
  double **dstatarray; /* array of size n_terms; the ith element in this
			  array is a pointer to an array of size
			  termarray[i].nstats                    */
} Model;

Model* ModelInitialize (char *fnames, char *sonames, double *inputs,
			int n_terms);

void ModelDestroy(Model *m);

/* A Model object contains information about an entire ERGM, including the
   total numbers of terms, parameters, and statistics along with a pointer
   to an array of ModelTerm structures.  */
typedef struct DegreeBoundstruct {
  int attrcount;
  int fBoundDegByAttr;
  int *attribs;
  int *maxout;
  int *minout;
  int *maxin;
  int *minin;
} DegreeBound;

DegreeBound* DegreeBoundInitialize(int *attribs, int *maxout, int *maxin, 
			   int *minout, int *minin, int condAllDegExact, 
			   int attriblength, Network *nwp);

void DegreeBoundDestroy(DegreeBound *bd);

int GetIndexForAttrValue(int value);
int ModelTermHamming (char *fnames, int n_terms);
int ModelTermFormation (char *fnames, int n_terms);
int ModelTermDissolve (char *fnames, int n_terms);

#endif

