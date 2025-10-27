/*  File inst/include/ergm_weighted_population.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2022 Statnet Commons
 */

#ifndef _ERGM_WEIGHTED_POPULATION_H_
#define _ERGM_WEIGHTED_POPULATION_H_

#include <R.h>

/* This is a data structure for weighted sampling.  It supports modifying weights 
   with relative efficiency when using the 'B' type; modifying weights when using 
   the 'W' type results in re-initializing the full data structure. */

typedef struct {
  char type;
  
  // binary stuff
  int height;
  double **weights;
  
  // Walker stuff
  int size;
  double *prob;
  int *alias;
  double *originalweights;
  double sum;
} WtPop;

static inline WtPop *WtPopInitialize(int size, double *weights, char type) {
  WtPop *wtp = R_Calloc(1, WtPop);
  
  if(size < 1) {
    error("cannot initialize weighted population of size less than 1");
  }

  for(int i = 0; i < size; i++) {
    if(weights[i] < 0) {
      error("cannot initialize weighted population with negative weights");
    }
  }
    
  if(type == 'B') {
    wtp->type = 'B';
    wtp->height = ceil(log2(size));

    wtp->weights = R_Calloc(wtp->height + 1, double *);
    for(int i = 0; i <= wtp->height; i++) {
      wtp->weights[i] = R_Calloc(pow(2,i), double);  
    }
    memcpy(wtp->weights[wtp->height], weights, size*sizeof(double));
    
    for(int i = wtp->height - 1; i >= 0; i--) {
      for(int j = pow(2,i) - 1; j >= 0; j--) {
        wtp->weights[i][j] = wtp->weights[i + 1][2*j] + wtp->weights[i + 1][2*j + 1];
      }
    }
    
    if(wtp->weights[0][0] == 0) {
      error("cannot initialize weighted population with zero total weight");
    }
  } else if(type == 'W') {
    wtp->type = 'W';
    wtp->size = size;      

    wtp->originalweights = R_Calloc(wtp->size, double);
    wtp->prob = R_Calloc(wtp->size, double);
    wtp->alias = R_Calloc(wtp->size, int);
    
    memcpy(wtp->originalweights, weights, wtp->size*sizeof(double));
    memcpy(wtp->prob, weights, wtp->size*sizeof(double));
    
    wtp->sum = 0;
    for(int i = 0; i < wtp->size; i++) {
      wtp->sum += wtp->prob[i];
    }

    if(wtp->sum == 0) {
      error("cannot initialize weighted population with zero total weight");
    }
        
    for(int i = 0; i < wtp->size; i++) {
      wtp->prob[i] = wtp->size*wtp->prob[i]/wtp->sum;
      wtp->alias[i] = -1;
    }
    
    // three passes to initialize;
    // underfulls and overfulls may initially occur in any order;
    // after the first pass, all underfulls precede all overfulls, 
    // and this is preserved during the second pass itself;
    // thus, after the second pass, all that can be left is roundoff 
    // error or categories that have prob exactly 1, and we set prob to 
    // 1 in either case (and also initialize alias for completeness)
    int i = 0;
    for(int pass = 0; pass < 2; pass++) {
      for(int j = 0; j < wtp->size; j++) {
        if(wtp->prob[j] < 1 && wtp->alias[j] < 0) {
          while(i < wtp->size && wtp->prob[i] <= 1) {
            i++;
          }
          if(i >= wtp->size) {
            break;
          }
          wtp->alias[j] = i;
          wtp->prob[i] -= 1 - wtp->prob[j];
        }
      }
    }
    
    for(int j = 0; j < wtp->size; j++) {
      if(wtp->alias[j] < 0) {
        wtp->alias[j] = j;
        wtp->prob[j] = 1;
      }
    }
  } else {
    error("unsupported weighted population type; options are 'B' for binary tree and 'W' for Walker");      
  }
  
  return wtp;
}

static inline void WtPopDestroy(WtPop *wtp) {
  if(wtp->type == 'B') {
    for(int i = 0; i <= wtp->height; i++) {
      R_Free(wtp->weights[i]);
    }
    R_Free(wtp->weights);
  } else {
    R_Free(wtp->prob);
    R_Free(wtp->alias);
    R_Free(wtp->originalweights);
  }
  
  R_Free(wtp);
}

static inline int WtPopGetRand(WtPop *wtp) {
  if(wtp->type == 'B') {
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
  } else {
    double s = unif_rand()*wtp->size;
    int i = (int) s;
    double y = s - i;
    if(y < wtp->prob[i]) {
      return i;  
    } else {
      return wtp->alias[i];
    }      
  }
}

static inline void WtPopSetWt(int position, double weight, WtPop *wtp) {
  if(wtp->type == 'B') {
    double change = weight - wtp->weights[wtp->height][position];
    int j = position;
    for(int i = wtp->height; i >= 0; i--) {
      wtp->weights[i][j] += change;
      j /= 2;
    }
  } else {
    wtp->originalweights[position] = weight;    
    WtPop *new_wtp = WtPopInitialize(wtp->size, wtp->originalweights, 'W');
    R_Free(wtp->prob);
    R_Free(wtp->alias);
    R_Free(wtp->originalweights);
    memcpy(wtp, new_wtp, sizeof(WtPop));
    R_Free(new_wtp);
  }
}

static inline double WtPopGetWt(int position, WtPop *wtp) {
  if(wtp->type == 'B') {
    return wtp->weights[wtp->height][position];
  } else {
    return wtp->originalweights[position];
  }
}

static inline double WtPopSumWts(WtPop *wtp) {
  if(wtp->type == 'B') {
    return wtp->weights[0][0];
  } else {
    return wtp->sum;
  }
}

#endif // _ERGM_WEIGHTED_POPULATION_H_
