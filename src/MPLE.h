#ifndef MPLE_H
#define MPLE_H

#include "edgetree.h"
#include "changestat.h"
#include "model.h"

void MPLE_wrapper (int *tails, int *heads, int *dnedges,
       int *maxpossibleedges,
		   int *dn, int *dflag, int *bipartite, int *nterms, 
		   char **funnames, char **sonames, double *inputs,  
		   int *responsevec, double *covmat,
		   int *weightsvector,
		   double * offset, double * compressedOffset,
		   int *maxNumDyadTypes, int *maxMPLEsamplesize, 
       int *compressflag);
void MpleInit_hash (int *responsevec, double *covmat, int *weightsvector,
		    double *offset, double *compressedOffset,
		    int maxNumDyadTypes, Edge maxMPLE, Network *nwp, Model *m);
void MpleInit_no_compress (int *responsevec, double *covmat, int *weightsvector,
		    double *offset, double *compressedOffset,
		    int maxNumDyadTypes, Edge maxMPLE, Network *nwp, Model *m);
int findCovMatRow(double *newRow,double *matrix, int rowLength, 
		  int numRows, int *responsevec, 
		  double * offset, double * compressedOffset, 
		  int curDyadNum); 
#endif
