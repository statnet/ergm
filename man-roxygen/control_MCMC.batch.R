#  File man-roxygen/control_MCMC.batch.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
#' @param MCMC.batch if not 0 or `NULL`, sample about this many
#'   networks per call to the lower-level code; this can be useful if
#'   `output=` is a function, where it can be used to limit the number
#'   of networks held in memory at any given time.
