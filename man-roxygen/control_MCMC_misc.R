#' @param seed Seed value (integer) for the random number generator.  See
#' \code{\link[base]{set.seed}}.
#' @param parallel Number of threads in which to run the sampling. Defaults to
#' 0 (no parallelism). See the entry on [parallel processing][ergm-parallel]
#' for details and troubleshooting.
#' @param parallel.type API to use for parallel processing. Supported values
#' are \code{"MPI"} and \code{"PSOCK"}. Defaults to using the \code{parallel}
#' package with PSOCK clusters. See \code{\link{ergm-parallel}}
#' @param parallel.version.check Logical: If TRUE, check that the version of
#' \code{\link[=ergm-package]{ergm}} running on the slave nodes is the same as
#' that running on the master node.
#' @param MCMC.packagenames Names of packages in which to look for change
#' statistic functions in addition to those autodetected. This argument should
#' not be needed outside of very strange setups.
