#ifndef _ERGM_OMP_H_
#define _ERGM_OMP_H_

#ifdef _OPENMP

int ergm_omp_terms;

#include <omp.h>
#define STRINGIFY(s) #s
#define ergm_PARALLEL_FOR _Pragma("omp parallel for if(ergm_omp_terms!=0) num_threads(ergm_omp_terms!=-1? ergm_omp_terms : omp_get_num_procs())")    
#define ergm_PARALLEL_FOR_LIMIT(lim) _Pragma(STRINGIFY(omp parallel for if(ergm_omp_terms!=0 && lim !=1) num_threads(MIN( lim ,ergm_omp_terms!=-1? ergm_omp_terms : omp_get_num_procs()))))
#else
#define ergm_PARALLEL_FOR
#define ergm_PARALLEL_FOR_LIMIT(lim)
#endif // _OMP

#endif // _ERGM_OMP_H_
