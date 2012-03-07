/*
 *  File ergm/src/plinfo.h
 *  Part of the statnet package, http://statnet.org
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) in
 *    http://statnet.org/attribution
 *
 *  Copyright 2012 the statnet development team
 */
#ifndef PLINFO_H
#define PLINFO_H

#include "edgetree.h"
#include "model.h"
#include "MPLE.h"

void plinfo_wrapper (int *tails, int *heads, int *dnedges,
		     int *dn, int *dflag, int *nterms, char **funnames,
		     char **sonames, double *inputs,  
		     double *responsevec, double *covmat,
		     int *start, int *end);
void plinfoInitialize (double *responsevec, double *covmat,
                       Vertex *start, Vertex *end,
                       Network *nwp, Model *m);

#endif
