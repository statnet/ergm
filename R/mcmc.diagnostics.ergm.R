#  File R/mcmc.diagnostics.ergm.R in package ergm, part of the Statnet suite of
#  packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
#=================================================================================
# This file contains the following 10 diagnostic tools and their helper functions
#      <mcmc.diagnostics>            <traceplot.ergm>
#      <mcmc.diagnostics.default>    <set.mfrow>
#      <mcmc.diagnostics.ergm>       <nvar.mcmc>
#      <is.mcmc.object>
#      <is.mcmc.list.object>
#      <plot.mcmc.ergm>              <varnames.mcmc>
#=================================================================================


#' Conduct MCMC diagnostics on a model fit
#'
#' This function prints diagnistic information and creates simple diagnostic
#' plots for MCMC sampled statistics produced from a fit.
#'
#' A pair of plots are produced for each statistic:a trace of the sampled
#' output statistic values on the left and density estimate for each variable
#' in the MCMC chain on the right.  Diagnostics printed to the console include
#' correlations and convergence diagnostics.
#'
#' \insertNoCite{RaLe95n}{ergm}
#'
#' @aliases mcmc.diagnostics.default
#' @param object A model fit object to be diagnosed.
#' @param \dots Additional arguments, to be passed to plotting functions.
#' @seealso [ergm()], \CRANpkg{network} package, \CRANpkg{coda} package,
#' [summary.ergm()]
#' @references \insertAllCited{}
#' @keywords models
#' @examples
#'
#' \dontrun{
#' #
#' data(florentine)
#' #
#' # test the mcmc.diagnostics function
#' #
#' gest <- ergm(flomarriage ~ edges + kstar(2))
#' summary(gest)
#'
#' #
#' # Plot the probabilities first
#' #
#' mcmc.diagnostics(gest)
#' #
#' # Use coda directly
#' #
#' library(coda)
#' #
#' plot(gest$sample, ask=FALSE)
#' #
#' # A full range of diagnostics is available
#' # using codamenu()
#' #
#' }
#'
#' @export mcmc.diagnostics
mcmc.diagnostics <- function(object, ...) {
  UseMethod("mcmc.diagnostics")
}

#' @noRd
#' @export
mcmc.diagnostics.default <- function(object, ...) {
  stop("An ergm object must be given as an argument ")
}

#' @describeIn mcmc.diagnostics
#'
#' @details For [ergm()] specifically, recent changes in the
#'   estimation algorithm mean that these plots can no longer be used
#'   to ensure that the mean statistics from the model match the
#'   observed network statistics. For that functionality, please use
#'   the GOF command: \code{gof(object, GOF=~model)}.
#'
#'   In fact, an [ergm()] output object contains the sample of
#'   statistics from the last MCMC run as element `$sample`. If
#'   missing data MLE is fit, the corresponding element is named
#'   `$sample.obs`. These are objects of [`mcmc`] and can be used
#'   directly in the \CRANpkg{coda} package to assess MCMC
#'   convergence.
#'
#'   More information can be found by looking at the documentation of
#'   [ergm()].
#'
#' @param center Logical: If `TRUE`, center the samples on the observed
#'   statistics.
#' @param esteq Logical: If `TRUE`, for statistics corresponding to
#'   curved ERGM terms, summarize the curved statistics by their
#'   negated estimating function values (evaluated at the MLE of any curved
#'   parameters) (i.e., \eqn{\eta'_{I}(\hat{\theta})\cdot (g_{I}(Y)-g_{I}(y))}
#'   for \eqn{I} being indices of the canonical parameters in
#'   question), rather than the canonical (sufficient) vectors of the
#'   curved statistics relative to the observed (\eqn{g_{I}(Y)-g_{I}(y)}).
#' @param vars.per.page Number of rows (one variable per row) per
#'   plotting page.  Ignored if \CRANpkg{latticeExtra} package is not
#'   installed.
#' @param which A character vector specifying which diagnostics to
#'   plot and/or print. Defaults to all of the below if meaningful:
#'   \describe{
#'
#' \item{`"plots"`}{Traceplots and density plots of sample values for all statistic or estimating function elements.}
#'
#' \item{`"texts"`}{Shorthand for the following text diagnostics.}
#'
#' \item{`"summary"`}{Summary of network statistic or estimating function elements as produced by [coda::summary.mcmc.list()].}
#'
#' \item{`"autocorrelation"`}{Autocorrelation of each of the network statistic or estimating function elements.}
#'
#' \item{`"crosscorrelation"`}{Cross-correlations between each pair of the network statistic or estimating function elements.}
#'
#' \item{`"burnin"`}{Burn-in diagnostics, in particular, the Geweke test.}
#'
#' } Partial matching is supported. (E.g., `which=c("auto","cross")`
#' will print autocorrelation and cross-correlations.)
#'
#' @param compact Numeric: For diagnostics that print variables in
#'   columns (e.g. correlations, hypothesis test p-values), try to
#'   abbreviate variable names to this many characters and round the
#'   numbers to `compact - 2` digits after the decimal point; 0 or
#'   `FALSE` for no abbreviation.
#'
#' @import coda
#' @export
mcmc.diagnostics.ergm <- function(object,
                                  center=TRUE,
                                  esteq=TRUE,
                                  vars.per.page=3,
                                  which=c("plots", "texts", "summary", "autocorrelation", "crosscorrelation", "burnin"),
                                  compact = FALSE, ...) {

  which <- match.arg(which, several.ok=TRUE)
  if("texts" %in% which) which <- c(which, "summary", "autocorrelation", "crosscorrelation", "burnin")

  # Coerce sample objects to mcmc.list. This allows all subsequent
  # operations to assume mcmc.list. The reason [["sample"]] is being
  # used here rather than $sample is because there is an unlikely
  # possibility that $sample doesn't exist but $sample.obs does.
  sm <- NVL3(object[["sample"]], as.mcmc.list(.))
  sm_thin <- attr(object[["sample"]], "extra_thin")
  sm.obs <- NVL3(object[["sample.obs"]], as.mcmc.list(.))
  sm_thin.obs <- attr(object[["sample.obs"]], "extra_thin")

  if(is.null(sm)) stop("MCMC was not run or MCMC sample was not stored.")

  if(!center){
    sm <- sweep.mcmc.list(sm, NVL(object$target.stats,object$nw.stats), "+")
    if(!is.null(sm.obs)) sm.obs <- sweep.mcmc.list(sm.obs, NVL(object$target.stats,object$nw.stats), "+")
  }else{
    # Then sm is already centered, *unless* there is missing data.  In
    # that case, center sm relative to sm.obs and center sm.obs to
    # 0. (The reason sm.obs is used as the reference is that it yields
    # the same result as if the whole network were observed.)
    if(!is.null(sm.obs)){
      sm.obs.mean <- colMeans.mcmc.list(sm.obs)
      sm <- sweep.mcmc.list(sm, sm.obs.mean, "-")
      sm.obs <- sweep.mcmc.list(sm.obs, sm.obs.mean, "-")
    }
  }

  if(esteq){
    if (!is.null(coef(object)) && !is.null(object$etamap)) {
      sm <- ergm.estfun(sm, theta=coef(object), model=object$etamap) %>% lapply.mcmc.list(`-`)
      if(!is.null(sm.obs)) sm.obs <- ergm.estfun(sm.obs, theta=coef(object), model=object$etamap) %>% lapply.mcmc.list(`-`)
    }
  }

  if("plots" %in% which){
    if(requireNamespace('latticeExtra', quietly=TRUE)){
      print(ergm_plot.mcmc.list(sm,main="Sample statistics",vars.per.page=vars.per.page,...))
      if(!is.null(sm.obs)) print(ergm_plot.mcmc.list(sm.obs,main="Constrained sample statistics",vars.per.page=vars.per.page,...))
    }else{
      plot(sm,...)
      if(!is.null(sm.obs)) plot(sm.obs,...)
    }
  }

  if("summary" %in% which){
    cat("Sample statistics summary:\n")
    print(summary(sm))

    if(!is.null(sm.obs)){
      cat("Constrained sample statistics summary:\n")
      print(summary(sm.obs))
    }
  }

  if(compact){
    varnames(sm) <- abbreviate(varnames(sm), compact)

    if(!is.null(sm.obs)){
      varnames(sm.obs) <- abbreviate(varnames(sm.obs), compact)
    }
  }else compact <- getOption("digits")+2

  if("summary" %in% which){
    # only show if we are using Hotelling termination criterion
    if(EVL(object$control$MCMLE.termination %in% c("Hotelling","precision","confidence"), TRUE)){
      # This can probably be improved.
      if(is.null(sm.obs)){
        cat("\nAre sample statistics significantly different from observed?\n")
        ds <- colMeans.mcmc.list(sm) - if(!center) NVL(object$target.stats, object$nw.stats) else 0
        sds <- apply(as.matrix(sm),2,sd)
        ns <- effectiveSize(sm)

        cv <-  cov(as.matrix(sm))

        z <- ds/sds*sqrt(ns)

      }else{
        cat("\nAre unconstrained sample statistics significantly different from constrained?\n")
        ds <- colMeans.mcmc.list(sm) - if(!center) colMeans.mcmc.list(sm.obs) else 0
        sds <- apply(as.matrix(sm),2,sd)
        sds.obs <- apply(as.matrix(sm.obs),2,sd)
        ns <- effectiveSize(sm)
        # It's OK constrained sample doesn't vary. (E.g, the extreme case
        # --- completely observed network --- is just one configuration of
        # statistics.)
        # Thus, the effective sample size for nonvarying is set to 1.
        ns.obs <- pmax(effectiveSize(sm.obs),1)

        cv <-  cov(as.matrix(sm))
        cv.obs <-  cov(as.matrix(sm.obs))

        z <- ds/sqrt(sds^2/ns+sds.obs^2/ns.obs)
      }
      p.z <- pnorm(abs(z),lower.tail=FALSE)*2

      overall.test <- approx.hotelling.diff.test(sm,sm.obs,if(is.null(sm.obs) && !center) object$target.stats else NULL)

      m <- rbind(c(ds,NA),c(z,overall.test$statistic),c(p.z,overall.test$p.value))
      rownames(m) <- c("diff.","test stat.","P-val.")
      colnames(m) <- c(varnames(sm), if(compact) "(Omni)" else "Omnibus (chi^2)")
      print(m, digits=compact-2)
    }
    # End simulated vs. observed test.
  }

  if("crosscorrelation" %in% which){
    cat("\nSample statistics cross-correlations:\n")
    print(crosscorr(sm), digits=compact-2)
    if(!is.null(sm.obs)){
      cat("Constrained sample statistics cross-correlations:\n")
      print(crosscorr(sm.obs), digits=compact-2)
    }

    cat("\nSample statistics auto-correlation:\n")
    for(chain in seq_along(sm)){
      ac<-autocorr(sm[[chain]],0:5)
      ac<-sapply(seq_len(ncol(sm[[chain]])),
                 function(i) ac[,i,i])
      colnames(ac)<-varnames(sm)
      cat("Chain", chain, "\n")
      print(ac, digits=compact-2)
    }
    if(!is.null(sm.obs)){
      cat("Constrained sample statistics auto-correlation:\n")
      for(chain in seq_along(sm.obs)){
        ac<-autocorr(sm.obs[[chain]],0:5)
        ac<-sapply(seq_len(ncol(sm.obs[[chain]])),
                   function(i) ac[,i,i])
        colnames(ac)<-varnames(sm.obs)
        cat("Chain", chain, "\n")
        print(ac, digits=compact-2)
      }
    }
  }

  if("burnin" %in% which){
    cat("\nSample statistics burn-in diagnostic (Geweke):\n")
    sm.gw<-geweke.diag(sm)
    sm.gws<-try(geweke.diag.mv(sm, split.mcmc.list=TRUE))
    if(!("try-error" %in% class(sm.gws))){
      for(chain in seq_along(sm.gw)){
        cat("Chain", chain, "\n")
        print(sm.gw[[chain]], digits=compact-2)
        cat("Individual P-values (lower = worse):\n")
        print(2*pnorm(abs(sm.gw[[chain]]$z),lower.tail=FALSE), digits=compact-2)
        cat("Joint P-value (lower = worse): ", format(sm.gws[[chain]]$p.value, digits=compact-2),"\n")
      }
    }
    if(!is.null(sm.obs)){
      cat("Sample statistics burn-in diagnostic (Geweke):\n")
      sm.obs.gw<-geweke.diag(sm.obs)
      sm.obs.gws<-try(geweke.diag.mv(sm.obs, split.mcmc.list=TRUE))
      if(!("try-error" %in% class(sm.obs.gws))){
        for(chain in seq_along(sm.obs.gw)){
          cat("Chain", chain, "\n")
          print(sm.obs.gw[[chain]], digits=compact-2)
          cat("P-values (lower = worse):\n")
          print(2*pnorm(abs(sm.obs.gw[[chain]]$z),lower.tail=FALSE), digits=compact-2)
          cat("Joint P-value (lower = worse): ", format(sm.gws[[chain]]$p.value, digits=compact-2),".\n")
        }
      }
    }
  }

  thin_note <- function(which, x, x_thin) if(NVL(x_thin, 1) != 1) writeLines(c("", strwrap(sprintf("Note: To save space, only one in every %d iterations of the%s MCMC sample used for estimation was stored for diagnostics. Sample size per chain was originally around %d with thinning interval %d.", x_thin, which, niter(x)*x_thin, thin(x)/x_thin), exdent = 2)))

  thin_note("", sm, sm_thin)
  thin_note(" constrained", sm.obs, sm_thin.obs)

  writeLines(c("", strwrap("Note: MCMC diagnostics shown here are from the last round of simulation, prior to computation of final parameter estimates. Because the final estimates are refinements of those used for this simulation run, these diagnostics may understate model performance. To directly assess the performance of the final model on in-model statistics, please use the GOF command: gof(ergmFitObject, GOF=~model).", exdent=2), ""))
}

#' Plot MCMC list using `lattice` package graphics
#'
#' @param x an [`mcmc.list`] object containing the mcmc diagnostic
#'   samples.
#' @param main character, main plot heading title.
#' @param vars.per.page Number of rows (one variable per row) per
#'   plotting page.  Ignored if \code{latticeExtra} package is not
#'   installed.
#' @param ... additional arguments, currently unused.
#' @note This is not a method at this time.
#'
#' @export ergm_plot.mcmc.list
ergm_plot.mcmc.list <- function(x, main=NULL, vars.per.page=3,...){
  if(!requireNamespace('lattice', quietly=TRUE) ||
     !requireNamespace('latticeExtra', quietly=TRUE))
    stop("ergm_plot.mcmc.list() requires ",sQuote('lattice')," and ",sQuote('latticeExtra')," packages.",call.=FALSE)

  # Workaround for coda::densityplot.mcmc.list(), which can't handle
  # duplicated variable names:
  varnames(x) <- make.unique(varnames(x))

  dp <- update(lattice::densityplot(x, panel=function(...){lattice::panel.densityplot(...);lattice::panel.abline(v=0)}, ...),xlab=NULL,ylab=NULL)
  tp <- update(lattice::xyplot(x, panel=function(...){lattice::panel.xyplot(...);lattice::panel.loess(...);lattice::panel.abline(0,0)}, ...),xlab=NULL,ylab=NULL)

  pages <- ceiling(nvar(x)/vars.per.page)
  # if the number of vars is less than vars.per.page, make adjustment
  if(nvar(x)<vars.per.page){
    vars.per.page<-nvar(x)
  }

  reordering <- c(rbind(seq_len(nvar(x)),nvar(x)+seq_len(nvar(x))))

  tpdp <- suppressWarnings(c(tp,dp))
  update(tpdp[reordering],layout=c(2,vars.per.page),as.table=TRUE,main=main)
}
