/*  File src/MPLE.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2019 Statnet Commons
 */
#ifndef MPLE_H
#define MPLE_H

#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_model.h"
#include "ergm_rlebdm.h"

void MPLE_wrapper(int *tails, int *heads, int *dnedges,
		  double *wl,
		  int *dn, int *dflag, int *bipartite, int *nterms, 
		  char **funnames, char **sonames, double *inputs,  
		  int *responsevec, double *covmat,
		  int *weightsvector,
		  int *maxDyads, int *maxDyadTypes);
void MpleInit_hash_wl_RLE(int *responsevec, double *covmat, int *weightsvector,
			  RLEBDM1D *wl, 
			  Edge maxDyads, Edge maxDyadTypes, Network *nwp, Model *m);
#endif
