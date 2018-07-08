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
  nnets <- length(unique(.peek_vattrv(nw, ".NetworkID")))

  if(is.numeric(subset)) subset <- unwhich(subset, nnets)
  subset <- rep(subset, length.out=nnets)
  remain <- subset

  stats <- stats.obs <- NULL
  
  while(any(remain)){
    message("Constructing GOF model.")
    pernet.m <- ergm_model(~N(GOF, subset=remain), nw=nw, response = object$response, ...,
                           term.options= modifyList(as.list(object$control$term.options), list(N.compact_stats=FALSE)))
    
    nstats <- nparam(pernet.m, canonical=TRUE)/sum(remain)
    
    message("Simulating unconstrained sample.")
    sim <- simulate(object, monitor=pernet.m, nsim=control$nsim, control=set.control.class("control.simulate.ergm",control), basis=nw, statsonly=TRUE, response = object$response, ...)
    sim <- sim[,ncol(sim)-sum(remain)*nstats+seq_len(sum(remain)*nstats),drop=FALSE]
  
    if(!is.null(object$constrained.obs)){
      message("Simulating constrained sample.")
      # FIXME: Simulations can be rerun only on the networks in the subset.
      sim.obs <- simulate(object, monitor=pernet.m, observational=TRUE, nsim=control$nsim, control=set.control.class("control.simulate.ergm",control), basis=nw, statsonly=TRUE, response = object$response, ...)
      sim.obs <- sim.obs[,ncol(sim.obs)-sum(remain)*nstats+seq_len(sum(remain)*nstats),drop=FALSE]
    }else{
      sim.obs <- matrix(summary(pernet.m, object$network, response = object$response, ...), nrow(sim), ncol(sim), byrow=TRUE)
    }

    message("Collating the simulations.")
    cn <- colnames(sim)[seq_len(nstats)] %>% sub(".*?:","", .)
    
    statarray <- array(c(sim), c(control$nsim, nstats, sum(remain)))
    dimnames(statarray) <- list(Iterations=NULL, Statistic=cn, Network=NULL)
    statarray.obs <- array(c(sim.obs), c(control$nsim, nstats, sum(remain)))
    dimnames(statarray.obs) <- list(Iterations=NULL, Statistic=cn, Network=NULL)

    if(is.null(stats)){
      stats <- statarray
      stats.obs <- statarray.obs
    }else{
      stats[,,remain[subset]] <- statarray
      stats.obs[,,remain[subset]] <- statarray.obs
    }

    # Calculate variances for each network and statistic.
    dv <- apply(statarray, 2:3, var) - apply(statarray.obs, 2:3, var)
    # If any statistic for the network has negative variance estimate, rerun it.
    remain[remain] <- apply(dv<=0, 2, any)
    if(any(remain)) message(sum(remain), " networks have bad simulations; rerunning.")
  }
  o <- structure(list(stats=stats, stats.obs=stats.obs, control=control), class="gofN")
}

#' @describeIn gofN A plotting method, making residuals vs. fitted and scale-location plots.
#' @export
plot.gofN <- function(x, ..., ask = dev.interactive()){
  if(ask){
    prev.ask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(prev.ask))
  }

  statarray <- x$stats
  statarray.obs <- x$stats.obs
  control <- x$control
  cn <- dimnames(statarray)[[2]]
  for(i in seq_along(cn)){
    s <- statarray[,i,]
    so <- statarray.obs[,i,]

    fitted <- colMeans(s)
    resid <- (colMeans(so)-colMeans(s))/sqrt(apply(s,2,var)-apply(so,2,var))

    plot(fitted, resid, ..., main = paste("Residuals vs. Fitted for", sQuote(cn[i])), xlab="Predicted value", ylab="Pearson residual",type="n")
    panel.smooth(fitted, resid, ...)
    abline(h=0, lty=3, col="gray")
    
    plot(fitted, sqrt(abs(resid)), ..., main = paste("Scale-location plot for", sQuote(cn[i])), xlab="Predicted value", ylab=expression(sqrt("|Pearson residual|")), type="n")
    panel.smooth(fitted, resid, ...)
    abline(h=0, lty=3, col="gray")
  }
}
