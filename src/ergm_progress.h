/*  File src/ergm_progress.h in package ergm, part of
 *  the Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#ifndef ERGM_PROGRESS_H
#define ERGM_PROGRESS_H

#include <R.h>
#include <Rinternals.h>

// Progress bar update frequency (update every N iterations)
#define ERGM_PROGRESS_UPDATE_FREQ 10
#define ERGM_PROGRESS_MPLE_UPDATE_FREQ 1000  // For large MPLE loops

// Progress bar functions that call R functions using cli package
// These are safe to call from C code

// Initialize a progress bar
// name: Progress bar label/message
// total: Total number of steps
// Returns: TRUE if successful, FALSE otherwise
static inline Rboolean ergm_progress_init(const char *name, unsigned int total) {
  SEXP call, result;
  int error = 0;
  
  // Protect all SEXP objects
  PROTECT(call = lang4(install(".ergm_progress_init"), 
                       mkString(name), 
                       ScalarInteger(total),
                       ScalarLogical(FALSE)));
  
  result = R_tryEval(call, R_GlobalEnv, &error);
  UNPROTECT(1);
  
  return error == 0;
}

// Update progress bar to current step
// current: Current step number (0-based or 1-based depending on usage)
static inline void ergm_progress_update(unsigned int current) {
  SEXP call;
  int error = 0;
  
  PROTECT(call = lang2(install(".ergm_progress_update"), 
                       ScalarInteger(current)));
  
  R_tryEval(call, R_GlobalEnv, &error);
  UNPROTECT(1);
}

// Complete and close the progress bar
static inline void ergm_progress_done(void) {
  SEXP call;
  int error = 0;
  
  PROTECT(call = lang1(install(".ergm_progress_done")));
  
  R_tryEval(call, R_GlobalEnv, &error);
  UNPROTECT(1);
}

#endif // ERGM_PROGRESS_H
