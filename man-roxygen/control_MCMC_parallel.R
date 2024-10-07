#  File man-roxygen/control_MCMC_parallel.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################
#' @param parallel Number of threads in which to run the sampling. Defaults to
#' 0 (no parallelism). See the entry on [parallel processing][ergm-parallel]
#' for details and troubleshooting.
#' @param parallel.type API to use for parallel processing. Supported values
#' are \code{"MPI"} and \code{"PSOCK"}. Defaults to using the \code{parallel}
#' package with PSOCK clusters. See \code{\link{ergm-parallel}}
#' @param parallel.version.check Logical: If TRUE, check that the version of
#' \pkg{ergm} running on the slave nodes is the same as
#' that running on the master node.
#' @param parallel.inherit.MT Logical: If TRUE, slave nodes and
#' processes inherit the [set.MT_terms()] setting.
