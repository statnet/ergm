#' @param MCMC.batch if not 0 or `NULL`, sample about this many
#'   networks per call to the lower-level code; this can be useful if
#'   `output=` is a function, where it can be used to limit the number
#'   of networks held in memory at any given time.
