#ifndef MPLE_H
#define MPLE_H

#include "edgeTree.h"
#include "basechangeStats.h"
#include "model.h"

void MPLE_wrapper (double *heads, double *tails, double *dnedges,
		   double *dn, int *dflag, double *bipartite, int *nterms, 
       char **funnames, char **sonames, double *inputs,  
		   int *responsevec, double *covmat,
		   int *weightsvector, double * offset, 
		   double * compressedOffset, int maxNumDyadTypes);
void MpleInitialize (int *responsevec, double *covmat,
	             int *weightsvector,
		     double * offset, double * compressedOffset,
		     int maxNumDyadTypes, Network *nwp, Model *m);
int rowsAreSame(double *rowA,double *rowB,int rowLength);
int findCovMatRow(double *newRow,double *matrix, int rowLength, 
		  int numRows, int *responsevec, 
		  double * offset, double * compressedOffset, 
		  int curDyadNum); 
void plinfo_wrapper (double *heads, double *tails, double *dnedges,
		     double *dn, int *dflag, int *nterms, char **funnames,
		     char **sonames, double *inputs,
		     double *responsevec, double *covmat,
		     int *start, int *end);
void plinfoInitialize (double *responsevec, double *covmat,
                       Vertex *start, Vertex *end,
                       Network *nwp, Model *m);
#endif
