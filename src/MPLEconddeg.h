/*
 *  File ergm/src/MPLEconddeg.h
 *  Part of the statnet package, http://statnet.org
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) in
 *    http://statnet.org/attribution
 *
 *  Copyright 2012 the statnet development team
 */
#ifndef MPLEconddeg_H
#define MPLEconddeg_H

#include "edgetree.h"
#include "changestat.h"
#include "MHproposal.h"
#include "model.h"

void MPLEconddeg_wrapper (int *tails, int *heads, int *dnedges,
                   int *dn, int *dflag, int *bipartite, 
                   int *nterms, char **funnames,
                   char **sonames, 
                   char **MHproposaltype, char **MHproposalpackage,
                   double *inputs, double *theta0, int *samplesize, 
                   double *sample, int *burnin, int *interval,  
                   int *newnetworktails, 
                   int *newnetworkheads, 
                   int *fVerbose, 
                   int *attribs, int *maxout, int *maxin, int *minout,
                   int *minin, int *condAllDegExact, int *attriblength, 
                   int *maxedges,
                   int *mtails, int *mheads, int *mdnedges);
void CondDegSampler (MHproposal *MHp,
		 double *theta, double *networkstatistics, 
		 int samplesize, int burnin, 
		 int interval, int fVerbose,
	       	 Network *nwp, Model *m);
void CondDegSample (MHproposal *MHp,
			 double *theta, double *statistics, 
			 int nsteps, int *staken,
			 int fVerbose,
			 Network *nwp, Model *m);
#endif
