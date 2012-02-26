#ifndef godfather_H 
#define godfather_H

#include <R.h>
#include "wtedgetree.h"
#include "model.h"

/* Function prototypes */
void godfather_wrapper (int *tails, int *heads, int *dnedges,
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

