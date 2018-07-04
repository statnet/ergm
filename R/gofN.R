#' Linear model diagnostics for multinetwork linear models
#'
#' @param object an [`ergm`] object.
#' @param GOF a one-sided [`ergm`] formula specifying network statistics whose goodness of fit to test.
#' @param subset argument for the [`N`][ergm-terms] term.
#' @param \dots additional arguments to `simulate.ergm()` and `summary.ergm_model()`. 
#' @param See [control.gof.ergm()].
#' @export
gofN <- function(object, GOF, subset=TRUE, control=control.gof.ergm(), ...){
  nw <- object$network

  pernet.m <- ergm_model(~N(GOF, subset=TRUE), nw=nw, response = object$response, ...,
                         term.options= modifyList(as.list(object$control$term.options), list(N.compact_stats=FALSE)))

  nnets <- length(unique(.peek_vattrv(nw, ".NetworkID")))
  nstats <- nparam(pernet.m, canonical=TRUE)/nnets

  sim <- simulate(object, monitor=pernet.m, nsim=control$nsim, control=set.control.class("control.simulate.ergm",control), statsonly=TRUE, response = object$response, ...)
  sim <- sim[,ncol(sim)-nnets*nstats+seq_len(nnets*nstats),drop=FALSE]
  
  if(!is.null(object$constrained.obs)){
    sim.obs <- simulate(object, monitor=pernet.m, constraints="obs", nsim=control$nsim, control=set.control.class("control.simulate.ergm",control), statsonly=TRUE, response = object$response, ...)
    sim.obs <- sim.obs[,ncol(sim.obs)-nnets*nstats+seq_len(nnets*nstats),drop=FALSE]
  }else{
    sim.obs <- matrix(summary(pernet.m, object$network, response = object$response, ...), nrow(sim), ncol(sim), byrow=TRUE)
  }

  cn <- colnames(sim)[seq_len(nstats)] %>% sub(".*?:","", .)
  
  statarray <- array(c(sim), c(control$nsim, nstats, nnets))
  dimnames(statarray) <- list(Iterations=NULL, Statistic=cn, Network=NULL)
  statarray.obs <- array(c(sim.obs), c(control$nsim, nstats, nnets))
  dimnames(statarray.obs) <- list(Iterations=NULL, Statistic=cn, Network=NULL)

  for(i in seq_along(cn)){
    s <- matrix(statarray[,i,],control$nsim,nnets)
    so <- matrix(statarray.obs[,i,],control$nsim,nnets)

    fitted <- colMeans(s)
    resid <- (colMeans(so)-colMeans(s))/sqrt(apply(s,2,var)-apply(so,2,var))

    plot(fitted, resid, ..., main = paste("Residuals vs. Fitted for", sQuote(cn[i])))
    abline(h=0, lty=3, col="gray")
    lines(lowess(fitted, resid), col="red")
  }
  invisible(list(stats=statarray, stats.obs=statarray.obs))
}
