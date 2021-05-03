#' @param
#'   MCMC.effectiveSize,MCMC.effectiveSize.damp,MCMC.effectiveSize.maxruns,MCMC.effectiveSize.burnin.pval,MCMC.effectiveSize.order.max
#'   Set `MCMC.effectiveSize` to a non-NULL value to adaptively
#'   determine the burn-in and the MCMC length needed to get the
#'   specified effective size; 50 is a reasonable value. In the
#'   adaptive MCMC mode, MCMC is run forward repeatedly
#'   (`MCMC.samplesize*MCMC.interval` steps, up to
#'   `MCMC.effectiveSize.maxruns` times) until the target effective
#'   sample size is reached or exceeded. After each run, the returned
#'   statistics are mapped to the estimating function scale, then a
#'   broken stick model is fit to each statistic to find the candidate
#'   burn-in. If its Geweke diagnostic produces a \eqn{p}-value higher
#'   than `MCMC.effectiveSize.burnin.pval`, it is accepted. The
#'   effective size of the post-burn-in sample is computed via Vats,
#'   Flegal, and Jones (2015), and compared to the target effective
#'   size. If it is not matched, the MCMC run is resumed, with the
#'   additional draws needed linearly extrapolated but weighted in
#'   favor of the baseline `MCMC.samplesize` by the weighting factor
#'   `MCMC.effectiveSize.damp` (higher = less damping). Lastly, if
#'   after an MCMC run, the number of samples equals or exceeds
#'   `2*MCMC.samplesize`, the chain will be thinned by 2 until it
#'   falls below that, while doubling
#'   `MCMC.interval`. `MCMC.effectiveSize.order.max` can be used to
#'   set the order of the AR model used to estimate the effective
#'   sample size and the variance for the Geweke diagnostic.
