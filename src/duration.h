#ifndef DURATION_H
#define DURATION_H

#include <R.h>
#include "edgeTree.h"
#include "MCMC.h"

/* Function prototypes */
void DurationMatrix (int *nedge, int *edge, int *ntimestep, int *nfem, 
      int *ntotal, int *nchange, int *change, int *ndissolve, int *dissolve,
      int *dmatrix);
void AddNewDurationRow (int *dmatrix, int row, int f, int m, int time, int offset);
void OverlapDurations (int *nedge, int *edge, int *ntimestep, int *nfem,
      int *ntotal, int *nchange, int *change, int *ndissolve, int *dissolve,
      int *maxoverlaps, int *omatrix);
void AddNewOverlapRow (int *omatrix, int row, int f1, int m1, 
int f2, int m2, int time1, int time2, int maxo);

void old_OverlapDurations (int *nedge, int *edge, int *ntimestep, int *nfem,
      int *ntotal, int *nchange, int *change, int *ndissolve, int *dissolve,
      int *maxoverlaps, int *omatrix);
void old_AddNewOverlapRow (int *omatrix, int row, int f1, int m1, 
int f2, int m2, int time, int maxo);



void godfather_wrapper (double *heads, double *tails, double *dnedges,
                   double *dn, int *dflag, double *bipartite, 
                   int *nterms, char **funnames,
                   char **sonames, 
                   double *totalntoggles, double *timestamps, 
                   double *toggleheads, double *toggletails,
                   double *inputs, 
                   double *changestats, 
                   int *newnetwork, 
                   int *fVerbose, 
                   double *maxedges);

      
#endif

