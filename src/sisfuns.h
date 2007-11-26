#ifndef SIS_H
#define SIS_H

#include "edgeTree.h"
#include "changestats.h"
#include "model.h"

void sisconj (int *colsums, int *conjseq, int nrow, int ncol, int i);
void sissort(int *rowsums, int *order, int nrow);
void sisknots (int *rowsums, int *conjseq, int *k, int *v, int nrow);
void sissamp(int *k, int *v, int *rowsums, int colsum, int nrow, int ncol, int *sampled, double *prob);
void sisufun(int *rowsums, int *sequence, int *ord, int nrow);
void removeknots(int *k, int *v, int n);
int get_max_range(int *list, int first, int last);
int get_min_range(int *list, int first, int last);
int sismin(int x, int y);
int sismax(int x, int y);
void SampleWithoutReplace(int start, int stop, double *weights, int nsample, int *sample, int nrow, double *prob);
void ComputeProb(double *weights, int nweights, int nsample, int weight_locat, double *prob);
void RecursiveProb(double *weights, int nweights, int nsample, double *prob);
//void sisstats(double *heads, double *tails, double *dnedges,
//		  double *dn, int *dflag, int *optionnum, char **funnames,
//		  char **sonames, double *inputs,  
//		  double *sample, 
//		  struct OptionInput *inp,
//		  int *fVerbose
//		);
//void sisdraw(int *rowsums, int *colsums, int *newmat, 
//	     int *nrow, int *ncol, 
//	     double *heads, double *tails,
//             double *dnedges, double *dn, int *dflag, int *optionnum, 
//	     char **funnames, char **sonames, 
//             double *inputs, 
//	     double *samplesize, double *graphstatistics,
//	     int *verb, double *prob, double *probvec,
//             double *params);
void sissim(int *rowsums, int *colsums, int *newmat,
	   int *nrow, int *ncol, 
	   double *samplesize, 
	   int *verb, double *prob, double *probvec);

#endif
