#ifndef INFECTION_H 
#define INFECTION_H

#include <R.h>
#include "edgeTree.h"
#include "MCMC.h"

/* Function prototypes */
void Prevalence (int *nnodes,
      int *nedge, int *edge, int *ntimestep, int *nfem,
      int *ntotal, int *nchange, int *change, int *ndissolve, int *dissolve,
      double *betarate, int *infected, int *prev);
      
#endif

