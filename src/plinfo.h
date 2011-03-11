#ifndef PLINFO_H
#define PLINFO_H

#include "edgetree.h"
#include "model.h"
#include "MPLE.h"

void plinfo_wrapper (int *tails, int *heads, int *dnedges,
         int *maxpossibleedges,
		     int *dn, int *dflag, int *nterms, char **funnames,
		     char **sonames, double *inputs,  
		     double *responsevec, double *covmat,
		     int *start, int *end);
void plinfoInitialize (double *responsevec, double *covmat,
                       Vertex *start, Vertex *end,
                       Network *nwp, Model *m);

#endif
