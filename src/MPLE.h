#ifndef MPLE_H
#define MPLE_H

#include "edgetree.h"
#include "changestat.h"
#include "model.h"

void MPLE_wrapper(int *tails, int *heads, int *dnedges,
		  int *wl, int *lel,
		  int *dn, int *dflag, int *bipartite, int *nterms, 
		  char **funnames, char **sonames, double *inputs,  
		  int *responsevec, double *covmat,
		  int *weightsvector,
		  int *maxDyads, int *maxDyadTypes);

void MpleInit_hash_bl(int *responsevec, double *covmat, int *weightsvector,
		      int *bl, 
		      Edge maxDyads, Edge maxDyadTypes, Network *nwp, Model *m);

void MpleInit_hash_wl(int *responsevec, double *covmat, int *weightsvector,
		      int *wl,
		      Edge maxDyadTypes, Network *nwp, Model *m);

#endif
