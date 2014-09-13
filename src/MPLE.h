/*  File src/MPLE.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2013 Statnet Commons
 */
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
