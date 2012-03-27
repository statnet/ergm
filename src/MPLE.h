#ifndef MPLE_H
#define MPLE_H

#include "edgetree.h"
#include "changestat.h"
#include "model.h"

void MPLE_wrapper(int *tails, int *heads, int *dnedges,
		  int *bltails, int *blheads, int *dblnedges,
		  int *dn, int *dflag, int *bipartite, int *nterms, 
		  char **funnames, char **sonames, double *inputs,  
		  int *responsevec, double *covmat,
		  int *weightsvector,
		  int *maxNumDyadTypes);

void MpleInit_hash(int *responsevec, double *covmat, int *weightsvector,
		   int *bltails, int *blheads, Edge blnedges, 
		   Edge maxNumDyadTypes, Network *nwp, Model *m);
#endif
