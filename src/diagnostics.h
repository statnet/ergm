#ifndef DIAGNOSTICS_H 
#define DIAGNOSTICS_H

#include <R.h>
#include "wtedgetree.h"
#include "model.h"

/* Function prototypes */
void DurationMatrix (int *nedge, int *edge, int *ntimestep,
      int *ntotal, int *nchange, int *change,
      int *dmatrix);
void AddNewDurationRow (int *dmatrix, int row, int h, int t, int time, int offset);
void OverlapDurations (int *nnodes, int *nedge, int *edge, int *ntimestep, int *nfem,
      int *ntotal, int *nchange, int *change,
      int *maxoverlaps, int *omatrix);
void AddNewOverlapRow (int *omatrix, int row, int f1, int m1, 
      int f2, int m2, int time1, int time2, int maxo);
void DegreeMixMatrix (int *nnodes,
      int *nedge, int *edge, int *ntimestep, int *nfem,
      int *ntotal, int *nchange, int *change,
      int *degmixmat);
void godfather_wrapper (int *tails, int *heads, int *dnedges,
      int *maxpossibleedges,
			int *dn, int *dflag, int *bipartite, 
			int *nterms, char **funnames,
			char **sonames, 
			int *totalntoggles, int *timestamps, 
			int *toggletails, int *toggleheads,
			int *dstart, int *dend,
			double *inputs, 
			double *changestats, 
			int *newnetworktails, 
			int *newnetworkheads, 
			int *accumulate, 
			int *fVerbose, 
			int *maxedges);
      
#endif

