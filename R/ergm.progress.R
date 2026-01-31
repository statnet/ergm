#  File R/ergm.progress.R in package ergm, part of
#  the Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons

# Progress bar helpers for C code
# These functions are called from C to provide progress updates

#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done

# Environment to store progress bar state
.progress_env <- new.env(parent = emptyenv())

#' Initialize a progress bar
#' @param name Name/message for the progress bar
#' @param total Total number of steps
#' @param clear Whether to clear the progress bar when done
#' @noRd
.ergm_progress_init <- function(name, total, clear = FALSE) {
  if (!is.null(.progress_env$id)) {
    # Clean up any existing progress bar
    try(cli::cli_progress_done(id = .progress_env$id), silent = TRUE)
  }
  
  .progress_env$id <- cli::cli_progress_bar(
    name,
    total = total,
    clear = clear,
    .auto_close = FALSE
  )
  
  invisible(.progress_env$id)
}

#' Update progress bar
#' @param current Current step number
#' @noRd
.ergm_progress_update <- function(current) {
  if (!is.null(.progress_env$id)) {
    cli::cli_progress_update(set = current, id = .progress_env$id)
  }
  invisible(NULL)
}

#' Complete and close progress bar
#' @noRd
.ergm_progress_done <- function() {
  if (!is.null(.progress_env$id)) {
    cli::cli_progress_done(id = .progress_env$id)
    .progress_env$id <- NULL
  }
  invisible(NULL)
}
