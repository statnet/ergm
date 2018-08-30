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

  message("Constructing simulation model(s).")
  sim_settings <- simulate(object, monitor=NULL, nsim=control$nsim, control=set.control.class("control.simulate.ergm",control), basis=nw, statsonly=TRUE, response = object$response, ..., do.sim=FALSE)
  if(!is.null(object$constrained.obs)) sim.obs_settings <- simulate(object, monitor=NULL, observational=TRUE, nsim=control$nsim, control=set.control.class("control.simulate.ergm",control), basis=nw, statsonly=TRUE, response = object$response, ..., do.sim=FALSE)

  prev.remain <- NULL
  while(any(remain)){
    message("Constructing GOF model.")
    if(NVL(prev.remain!=remain, TRUE))
      pernet.m <- ergm_model(~N(GOF, subset=remain), nw=nw, response = object$response, ...,
                             term.options= modifyList(as.list(object$control$term.options), list(N.compact_stats=FALSE)))
    prev.remain <- remain
    
    nstats <- nparam(pernet.m, canonical=TRUE)/sum(remain)
    
    message("Simulating unconstrained sample.")
    sim <- do.call(simulate, modifyList(sim_settings, list(monitor=pernet.m)))
    sim <- sim[,ncol(sim)-sum(remain)*nstats+seq_len(sum(remain)*nstats),drop=FALSE]
  
    if(!is.null(object$constrained.obs)){
      message("Simulating constrained sample.")
      # FIXME: Simulations can be rerun only on the networks in the subset.
      sim.obs <- do.call(simulate, modifyList(sim.obs_settings, list(monitor=pernet.m)))
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
    if(any(remain)) message(sum(remain), " networks (", paste(which(remain),collapse=", "), ") have bad simulations; rerunning.")
  }
  o <- structure(list(nw=nw, subset=subset, stats=stats, stats.obs=stats.obs, control=control), class="gofN")
}

#' @describeIn gofN A plotting method, making residuals and scale-location plots.
#' 
#' @param against vector of values, network attribute, or a formula whose RHS gives an expression in terms of network attributes to plot against; if `NULL` (default), plots against fitted values.
#' @param col,pch,cex vector of values (wrapped in [I()]), network attribute, or a formula whose RHS gives an expression in terms of network attributes to plot against.
#' @param which which to plot (`1` for residuals plot, `2` for \eqn{\sqrt{|R_i|}} scale plot).
#' 
#' @export
plot.gofN <- function(x, against=NULL, which=1:2, col=1, pch=1, cex=1, ..., ask = dev.interactive()){
  if(ask){
    prev.ask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(prev.ask))
  }

  nattrs <- get_multinet_nattr_tibble(x$nw)[x$subset,]
  
  againstname <- switch(class(against),
                        character = against,
                        formula = despace(deparse(against[[length(against)]])),
                        `NULL` = "Predicted value",
                        despace(deparse(substitute(against))))
  againstval <- switch(class(against),
                       character = nattrs[[against]],
                       formula = eval(against[[length(against)]], envir = nattrs, enclos = environment(against)),
                       against)

  for(gpar in c("col", "pch", "cex")){
    a <- get(gpar)
    a <- switch(class(a),
                AsIs = a,
                character = nattrs[[a]],
                formula = eval(a[[length(a)]], envir = nattrs, enclos = environment(a)),
                a)
    assign(gpar, a)
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

    if(1L %in% which){
      plot(NVL(againstval,fitted), resid, col=col, pch=pch, cex=cex,..., main = paste("Residuals vs. Fitted for", sQuote(cn[i])), xlab=againstname, ylab="Pearson residual",type="n")
      panel.smooth(NVL(againstval,fitted), resid, col=col, pch=pch, cex=cex, ...)
      abline(h=0, lty=3, col="gray")
    }
    
    if(2L %in% which){
      plot(NVL(againstval,fitted), sqrt(abs(resid)), col=col, pch=pch, cex=cex,..., main = paste("Scale-location plot for", sQuote(cn[i])), xlab=againstname, ylab=expression(sqrt("|Pearson residual|")), type="n")
      panel.smooth(NVL(againstval,fitted), sqrt(abs(resid)), col=col, pch=pch, cex=cex, ...)
      abline(h=0, lty=3, col="gray")
    }
  }
}
