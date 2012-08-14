#ifndef MPLE_H
#define MPLE_H

#include "edgetree.h"
#include "changestat.h"
#include "model.h"

void MPLE_wrapper(int *tails, int *heads, int *dnedges,
		  int *wl, int *ltails, int *lheads, int *dlnedges,
		  int *dn, int *dflag, int *bipartite, int *nterms, 
		  char **funnames, char **sonames, double *inputs,  
		  int *responsevec, double *covmat,
		  int *weightsvector,
		  int *maxNumDyadTypes);

void MpleInit_hash_bl(int *responsevec, double *covmat, int *weightsvector,
		      int *bltails, int *blheads, Edge blnedges, 
		      Edge maxNumDyadTypes, Network *nwp, Model *m);

void MpleInit_hash_wl(int *responsevec, double *covmat, int *weightsvector,
		      int *wltails, int *wlheads, Edge wlnedges, 
		      Edge maxNumDyadTypes, Network *nwp, Model *m);

#endif
