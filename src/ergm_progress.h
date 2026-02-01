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
#include <cli/progress.h>

// Progress bar update frequency (update every N iterations)
#define ERGM_PROGRESS_UPDATE_FREQ 10
#define ERGM_PROGRESS_MPLE_UPDATE_FREQ 1000  // For large MPLE loops

// Progress bar structure to track state
typedef struct {
  SEXP id;  // Progress bar ID from cli
  unsigned int total;
  Rboolean active;
} ergm_progress_bar;

// Initialize a progress bar using cli C API
// name: Progress bar label/message
// total: Total number of steps
// Returns: Progress bar structure
static inline ergm_progress_bar ergm_progress_init(const char *name, unsigned int total) {
  ergm_progress_bar pb = {R_NilValue, total, FALSE};
  
  if (CLI_SHOULD_TICK) {
    // Create progress bar using cli C API
    pb.id = PROTECT(cli_progress_bar(total, NULL));
    pb.active = TRUE;
  }
  
  return pb;
}

// Update progress bar to current step
// pb: Progress bar structure
// current: Current step number
static inline void ergm_progress_update(ergm_progress_bar *pb, unsigned int current) {
  if (pb->active && !Rf_isNull(pb->id)) {
    cli_progress_set(pb->id, current);
  }
}

// Complete and close the progress bar
// pb: Progress bar structure
static inline void ergm_progress_done(ergm_progress_bar *pb) {
  if (pb->active && !Rf_isNull(pb->id)) {
    cli_progress_done(pb->id);
    UNPROTECT(1);  // Unprotect the id that was protected in init
    pb->active = FALSE;
  }
}

#endif // ERGM_PROGRESS_H
