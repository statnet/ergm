/*  File inst/include/ergm_Rutil.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2023 Statnet Commons
 */
#ifndef _ERGM_RUTIL_H_
#define _ERGM_RUTIL_H_
#include<Rinternals.h>
#include<string.h>

/*
  This function is based on
  https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Handling-lists

  I'm putting it here pending its terms of use being clarified in https://bugs.r-project.org/bugzilla/show_bug.cgi?id=17664 .
*/
static inline SEXP getListElement(SEXP list, const char *str){
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);

  for (unsigned int i = 0; i < length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  return elmt;
}

static inline SEXP setListElement(SEXP list, const char *str, SEXP value){
  SEXP names = getAttrib(list, R_NamesSymbol);

  for (unsigned int i = 0; i < length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      SET_VECTOR_ELT(list, i, value);
      return value;
    }
  error("List does not have element '%s' to set.", str);
  return R_NilValue;
}

static inline SEXP mkRStrVec(const char **x){
  unsigned int l=0;
  while(x[l]) l++;
  SEXP o = PROTECT(allocVector(STRSXP, l));
  for(unsigned int i=0; i<l; i++) SET_STRING_ELT(o, i, mkChar(x[i]));
  UNPROTECT(1);
  return o;
}

#define TOINTSXP(x) x = PROTECT(coerceVector(x, INTSXP))
#define TOREALSXP(x) x = PROTECT(coerceVector(x, REALSXP))
#define FIRSTCHAR(x) CHAR(STRING_ELT(x, 0))

// Safely test if x is a NULL C pointer or a NULL R pointer.
#define isNULL(x) ((x)==NULL || (x)==R_NilValue)

// An alias for R_alloc that behaves more like Calloc(); uses the
// following helper function:
static inline void *R_calloc_helper(size_t n, size_t size){
  char *tmp = R_alloc(n, size);
  memset(tmp, 0, n*size);
  return (void*) tmp;
}
#define R_calloc(n, type) ((type*) R_calloc_helper((n), sizeof(type)))

#define R_CheckUserInterruptEvery(freq, iter) if((iter)%(freq) == 0) R_CheckUserInterrupt();

#endif 
