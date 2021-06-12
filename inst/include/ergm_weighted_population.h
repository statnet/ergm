/*  File inst/include/ergm_weighted_population.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2021 Statnet Commons
 */
#ifndef _ERGM_WEIGHTED_POPULATION_H_
#define _ERGM_WEIGHTED_POPULATION_H_

#include <R.h>

/* This is a data structure for weighted sampling from a variable
   population.  For now, these are static-inlined for simplicity, 
   but we may want to move some of the less frequently used functions
   to a C file and save RAM. */

typedef struct{
  double **weights;
  int height;
} WtPop;

static inline WtPop *WtPopInitialize(int size, double *weights){
  WtPop *wtp = Calloc(1, WtPop);

  wtp->height = ceil(log2(size));
  wtp->weights = Calloc(wtp->height + 1, double *);
  for(int i = 0; i <= wtp->height; i++) {
    wtp->weights[i] = Calloc(pow(2,i), double);  
  }
  memcpy(wtp->weights[wtp->height], weights, size*sizeof(double));
  
  for(int i = wtp->height - 1; i >= 0; i--) {
    for(int j = pow(2,i) - 1; j >= 0; j--) {
      wtp->weights[i][j] = wtp->weights[i + 1][2*j] + wtp->weights[i + 1][2*j + 1];
    }
  }
  
  return wtp;
}

static inline void WtPopDestroy(WtPop *wtp){
  for(int i = 0; i <= wtp->height; i++) {
    Free(wtp->weights[i]);
  }
  Free(wtp->weights);
  Free(wtp);
}

static inline int WtPopGetRand(WtPop *wtp){
  double s = unif_rand()*wtp->weights[0][0];
  int j = 0;
  for(int i = 1; i <= wtp->height; i++) {
    if(s > wtp->weights[i][2*j]) {
      s -= wtp->weights[i][2*j];
      j = 2*j + 1;
    } else {
      j = 2*j;
    }
  }
  
  return j;
}

static inline void WtPopSetWt(int position, double weight, WtPop *wtp){
  double change = weight - wtp->weights[wtp->height][position];
  int j = position;
  for(int i = wtp->height; i >= 0; i--) {
    wtp->weights[i][j] += change;
    j /= 2;
  }
}

static inline double WtPopGetWt(int position, WtPop *wtp){
  return wtp->weights[wtp->height][position];
}

static inline double WtPopSumWts(WtPop *wtp) {
  return wtp->weights[0][0];
}

#endif // _ERGM_WEIGHTED_POPULATION_H_
