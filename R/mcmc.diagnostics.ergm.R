#  File R/mcmc.diagnostics.ergm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################
#=================================================================================
# This file contains the following 10 diagnostic tools and their helper functions
#      <mcmc.diagnostics>            <traceplot.ergm>
#      <mcmc.diagnostics.default>    <set.mfrow>
#      <mcmc.diagnostics.ergm>       <nvar.mcmc>
#      <is.mcmc.object>
#      <is.mcmc.list.object>
#      <plot.mcmc.ergm>              <varnames.mcmc>
#=================================================================================



#########################################################################
# The <mcmc.diagnostics.X> functions create diagnostic plots for the
# MCMC sampled statistics of the ergm X and prints the Raftery-Lewis
# diagnostics, indicating whether they are sufficient or not; if X is not
# an ergm, execution will halt
#
# --PARAMTERS--
#   object : an ergm object, that has an MCMC established stats matrix
#   sample : the name of the component in 'object' to base the diagnosis
#            on; recognized strings are "observed", "sample", and
#            "thetasample"; default="sample"
#   smooth : whether to draw a smooth line through the trace plots;
#            default=TRUE
#   maxplot: the maximum number of statistics to plot; default=1000
#   verbose: whether to print out additional information about the
#            MCMC runs including lag correlations; default=TRUE
#   center : whether the samples should be centered on the observed
#            statistics; default=TRUE
#   ...    : addtional parameters that are passed to <plot.mcmc.ergm>
#   main, ylab, xlab: have their usual par-like meanings
#
#
# --IGNORED PARAMETERS--
#   r      : what percentile of the distribution to estimate; this is
#            ignored: default=.0125
#   digits : the number of digits to print; default=6
#
# --RETURNED--
#   raft: a list containing
#    params          : ?? 
#    resmatrix       : ??
#    degeneracy.value: the degeneracy.value of 'object', as computed by
#                      <ergm.degeneracy>
#    degeneracy.type : the degeneracy.type of 'object', as computed by
#                      <ergm.compute.degeneracy> 
#    simvals         :
#
##########################################################################



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
#' @aliases mcmc.diagnostics.default
#' @param object A model fit object to be diagnosed.
#' @param \dots Additional arguments, to be passed to plotting functions.
#' @seealso \code{\link{ergm}}, \code{network} package, \code{coda} package,
#' \code{\link{summary.ergm}}
#' @references % Warnes, G.W. (2000).  Multi-Chain and Parallel Algorithms for
#' Markov % Chain Monte Carlo. Dissertation, Department of Biostatistics, %
#' University of Washington, % Raftery, A.E. and Lewis, S.M. (1992).  One long
#' run with diagnostics: Implementation strategies for Markov chain Monte
#' Carlo. Statistical Science, 7, 493-497.
#' 
#' Raftery, A.E. and Lewis, S.M. (1995).  The number of iterations, convergence
#' diagnostics and generic Metropolis algorithms.  In Practical Markov Chain
#' Monte Carlo (W.R. Gilks, D.J. Spiegelhalter and S. Richardson, eds.).
#' London, U.K.: Chapman and Hall.
#' 
#' This function is based on the \code{coda} package It is based on the the R
#' function \code{raftery.diag} in \code{coda}.  \code{raftery.diag}, in turn,
#' is based on the FORTRAN program \code{gibbsit} written by Steven Lewis which
#' is available from the Statlib archive.
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
#'   In fact, an ergm output \code{object} contains the matrix of
#'   statistics from the MCMC run as component \code{$sample}.  This
#'   matrix is actually an object of class \code{mcmc} and can be used
#'   directly in the \code{coda} package to assess MCMC
#'   convergence. \emph{Hence all MCMC diagnostic methods available in
#'   \code{coda} are available directly.} See the examples and
#'   \url{https://www.mrc-bsu.cam.ac.uk/software/bugs/the-bugs-project-winbugs/coda-readme/}.
#' 
#'   More information can be found by looking at the documentation of
#'   \code{\link{ergm}}.
#' 
#' @param center Logical: If TRUE, center the samples on the
#'   observed statistics.
#' @param esteq Logical: If TRUE, for statistics corresponding to
#'   curved ERGM terms, summarize the curved statistics by their
#'   estimating equation values (evaluated at the MLE of any curved
#'   parameters) (i.e., \eqn{\eta'_{I}(\hat{\theta})\cdot g_{I}(y)}
#'   for \eqn{I} being indices of the canonical parameters in
#'   question), rather than the canonical (sufficient) vectors of the
#'   curved statistics (\eqn{g_{I}(y)}).
#' @param vars.per.page Number of rows (one variable per row) per
#'   plotting page.  Ignored if \code{latticeExtra} package is not
#'   installed.
#' @return \code{\link{mcmc.diagnostics.ergm}} returns some degeneracy
#' information, if it is included in the original object.  The function is
#' mainly used for its side effect, which is to produce plots and summary
#' output based on those plots.
#' @import coda
#' @export
mcmc.diagnostics.ergm <- function(object,
                                  center=TRUE,
                                  esteq=TRUE,
                                  vars.per.page=3,...) {
#
  if(!is.null(object$degeneracy.value) && !is.na(object$degeneracy.value)){
   degeneracy.value <- object$degeneracy.value
   degeneracy.type <- object$degeneracy.type
  }else{
   degeneracy.value <- NULL
   degeneracy.type <- NULL
  }

  # Coerce sample objects to mcmc.list. This allows all subsequent
  # operations to assume mcmc.list. The reason [["sample"]] is being
  # used here rather than $sample is because there is an unlikely
  # possibility that $sample doesn't exist but $sample.obs does.
  sm <- NVL3(object[["sample"]], as.mcmc.list(.))
  sm.obs <- NVL3(object[["sample.obs"]], as.mcmc.list(.))

  if(is.null(sm)) stop("MCMC was not run or MCMC sample was not stored.")

  if(!center){
    sm <- sweep.mcmc.list(sm, object$target.stats, "+")
    if(!is.null(sm.obs)) sm.obs <- sweep.mcmc.list(sm.obs, object$target.stats, "+")
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
    if (!is.null(object$coef) && !is.null(object$etamap)) {
      sm <- ergm.estfun(sm, theta=coef(object), model=object$etamap)
      if(!is.null(sm.obs)) sm.obs <- ergm.estfun(sm.obs, theta=coef(object), model=object$etamap)
    }
  }

  cat("Sample statistics summary:\n")
  print(summary(sm))
  if(!is.null(sm.obs)){
    cat("Constrained sample statistics summary:\n")
    print(summary(sm.obs))
  }
  
  # only show if we are using Hotelling termination criterion
  if (identical(object$control$MCMLE.termination, "Hotelling")) {
    # This can probably be improved.
    if(is.null(sm.obs)){
      cat("\nAre sample statistics significantly different from observed?\n")
      ds <- colMeans.mcmc.list(sm) - if(!center) object$target.stats else 0
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
    colnames(m) <- c(varnames(sm),"Overall (Chi^2)")
    print(m)
  }
  # End simulated vs. observed test.
  
  cat("\nSample statistics cross-correlations:\n")
  print(crosscorr(sm))
  if(!is.null(sm.obs)){
    cat("Constrained sample statistics cross-correlations:\n")
    print(crosscorr(sm.obs))
  }

  cat("\nSample statistics auto-correlation:\n")
  for(chain in seq_along(sm)){
    ac<-autocorr(sm[[chain]],0:5)
    ac<-sapply(seq_len(ncol(sm[[chain]])),
               function(i) ac[,i,i])
    colnames(ac)<-varnames(sm)
    cat("Chain", chain, "\n")
    print(ac)
  }
  if(!is.null(sm.obs)){
    cat("Constrained sample statistics auto-correlation:\n")
      for(chain in seq_along(sm.obs)){
        ac<-autocorr(sm.obs[[chain]],0:5)
        ac<-sapply(seq_len(ncol(sm.obs[[chain]])),
                   function(i) ac[,i,i])
        colnames(ac)<-varnames(sm.obs)
        cat("Chain", chain, "\n")
        print(ac)
      }
  }

  cat("\nSample statistics burn-in diagnostic (Geweke):\n")
  sm.gw<-geweke.diag(sm)
  sm.gws<-try(geweke.diag.mv(sm, split.mcmc.list=TRUE))
  if(!("try-error" %in% class(sm.gws))){
  for(chain in seq_along(sm.gw)){
    cat("Chain", chain, "\n")
    print(sm.gw[[chain]])
    cat("Individual P-values (lower = worse):\n")
    print(2*pnorm(abs(sm.gw[[chain]]$z),lower.tail=FALSE))
    cat("Joint P-value (lower = worse): ", sm.gws[[chain]]$p.value,".\n")
  }
  }
  if(!is.null(sm.obs)){
    cat("Sample statistics burn-in diagnostic (Geweke):\n")
    sm.obs.gw<-geweke.diag(sm.obs)
    sm.obs.gws<-try(geweke.diag.mv(sm.obs, split.mcmc.list=TRUE))
    if(!("try-error" %in% class(sm.obs.gws))){
    for(chain in seq_along(sm.obs.gw)){
      cat("Chain", chain, "\n")
      print(sm.obs.gw[[chain]])
      cat("P-values (lower = worse):\n")
      print(2*pnorm(abs(sm.obs.gw[[chain]]$z),lower.tail=FALSE))
      cat("Joint P-value (lower = worse): ", sm.gws[[chain]]$p.value,".\n")
    }
   }
  }
  
  if(requireNamespace('latticeExtra', quietly=TRUE)){  
    print(ergm_plot.mcmc.list(sm,main="Sample statistics",vars.per.page=vars.per.page,...))
    if(!is.null(sm.obs)) print(ergm_plot.mcmc.list(sm.obs,main="Constrained sample statistics",vars.per.page=vars.per.page,...))
  }else{
    plot(sm,...)
    if(!is.null(sm.obs)) plot(sm.obs,...)
  }
  
  cat("\nMCMC diagnostics shown here are from the last round of simulation, prior to computation of final parameter estimates. Because the final estimates are refinements of those used for this simulation run, these diagnostics may understate model performance. To directly assess the performance of the final model on in-model statistics, please use the GOF command: gof(ergmFitObject, GOF=~model).\n")

  invisible(list(degeneracy.value=degeneracy.value,
                 degeneracy.type=degeneracy.type))
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
  if(!requireNamespace('lattice', quietly=TRUE, warn.conflicts=FALSE) ||
     !requireNamespace('latticeExtra', quietly=TRUE, warn.conflicts=FALSE))
    stop("ergm_plot.mcmc.list() requires ",sQuote('lattice')," and ",sQuote('latticeExtra')," packages.",call.=FALSE)
  
  dp <- update(lattice::densityplot(x, panel=function(...){lattice::panel.densityplot(...);lattice::panel.abline(v=0)}),xlab=NULL,ylab=NULL)
  tp <- update(lattice::xyplot(x, panel=function(...){lattice::panel.xyplot(...);lattice::panel.loess(...);lattice::panel.abline(0,0)}),xlab=NULL,ylab=NULL)
  
  pages <- ceiling(nvar(x)/vars.per.page)
  # if the number of vars is less than vars.per.page, make adjustment
  if(nvar(x)<vars.per.page){
    vars.per.page<-nvar(x)
  }
  
  reordering <- c(rbind(seq_len(nvar(x)),nvar(x)+seq_len(nvar(x))))
  
  tpdp <- suppressWarnings(c(tp,dp))
  update(tpdp[reordering],layout=c(2,vars.per.page),as.table=TRUE,main=main)
}
