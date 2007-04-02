#ifndef DURATION_H
#define DURATION_H

#include <R.h>
#include "edgeTree.h"


/* Function prototypes */
void DurationMatrix (int *nedge, int *edge, int *ntimestep, int *nfem, 
      int *ntotal, int *nchange, int *change, int *ndissolve, int *dissolve,
      int *dmatrix);
void AddNewDurationRow (int *dmatrix, int row, int f, int m, int time, int offset);
void OverlapDurations (int *nedge, int *edge, int *ntimestep, int *nfem,
      int *ntotal, int *nchange, int *change, int *ndissolve, int *dissolve,
      int *maxoverlaps, int *omatrix);
void AddNewOverlapRow (int *omatrix, int row, int f1, int m1, 
      int f2, int m2, int time, int maxo);

#endif

