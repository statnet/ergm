#  File man-roxygen/control_MCMC_effectiveSize.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
#' @param
#'   MCMC.effectiveSize,MCMC.effectiveSize.damp,MCMC.effectiveSize.maxruns,MCMC.effectiveSize.burnin.pval,MCMC.effectiveSize.burnin.min,MCMC.effectiveSize.burnin.max,MCMC.effectiveSize.burnin.nmin,MCMC.effectiveSize.burnin.nmax,MCMC.effectiveSize.burnin.PC,MCMC.effectiveSize.burnin.scl,MCMC.effectiveSize.order.max
#'   Set `MCMC.effectiveSize` to a non-NULL value to adaptively
#'   determine the burn-in and the MCMC length needed to get the
#'   specified effective size; 50 is a reasonable value. In the
#'   adaptive MCMC mode, MCMC is run forward repeatedly
#'   (`MCMC.samplesize*MCMC.interval` steps, up to
#'   `MCMC.effectiveSize.maxruns` times) until the target effective
#'   sample size is reached or exceeded.
#'
#'   After each run, the returned statistics are mapped to the
#'   estimating function scale, then an exponential decay model is fit
#'   to the scaled statistics to find that burn-in which would reduce
#'   the difference between the initial values of statistics and their
#'   equilibrium values by a factor of `MCMC.effectiveSize.burnin.scl`
#'   of what it initially was, bounded by `MCMC.effectiveSize.min` and
#'   `MCMC.effectiveSize.max` as proportions of sample size. If the
#'   best-fitting decay exceeds `MCMC.effectiveSize.max`, the
#'   exponential model is considered to be unsuitable and
#'   `MCMC.effectiveSize.min` is used.
#'
#'   A Geweke diagnostic is then run, after thinning the sample to
#'   `MCMC.effectiveSize.burnin.nmax`. If this Geweke diagnostic
#'   produces a \eqn{p}-value higher than
#'   `MCMC.effectiveSize.burnin.pval`, it is accepted.
#'
#'   If `MCMC.effectiveSize.burnin.PC>0`, instead of using the full
#'   sample for burn-in estimation, at most this many principal
#'   components are used instead.
#'
#'   The effective
#'   size of the post-burn-in sample is computed via Vats, Flegal, and
#'   Jones (2015), and compared to the target effective size. If it is
#'   not matched, the MCMC run is resumed, with the additional draws
#'   needed linearly extrapolated but weighted in favor of the
#'   baseline `MCMC.samplesize` by the weighting factor
#'   `MCMC.effectiveSize.damp` (higher = less damping). Lastly, if
#'   after an MCMC run, the number of samples equals or exceeds
#'   `2*MCMC.samplesize`, the chain will be thinned by 2 until it
#'   falls below that, while doubling
#'   `MCMC.interval`. `MCMC.effectiveSize.order.max` can be used to
#'   set the order of the AR model used to estimate the effective
#'   sample size and the variance for the Geweke diagnostic.
