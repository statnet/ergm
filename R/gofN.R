.mean_var <- function(x, s){
  ng <- length(x)/s
  split(x, rep(1:ng,each=s)) %>%
    map(var) %>%
    reduce(`+`) %>%
    `/`(.,ng)
}

#' Linear model diagnostics for multinetwork linear models
#'
#' @param object an [`ergm`] object.
#' @param GOF a one-sided [`ergm`] formula specifying network
#'   statistics whose goodness of fit to test.
#' @param subset argument for the [`N`][ergm-terms] term.
#' @param \dots additional arguments to [simulate.ergm()] and
#'   [summary.ergm_model()].
#' @param control See [control.gof.ergm()].
#' @param obs.twostage Either `FALSE` or an integer. If not `FALSE`,
#'   the constrained sample variance is estimated by simulating
#'   conditional on the observed networks. However, a more accurate
#'   estimate can be obtained via a two-stage process: \enumerate{
#' 
#' \item Sample networks from the model without the observational
#' constraint.
#'
#' \item Conditional on each of those networks, sample with the
#' observational constraint, estimating the variance within each
#' sample and then averaging over the first-stage sample.
#'
#' }
#'
#' Then, `obs.twostage` specifies the number of unconstrained networks
#' to simulate from, which should divide the [control.gof.ergm()]'s
#' `nsim` argument evenly.
#'
#' @export
gofN <- function(object, GOF, subset=TRUE, control=control.gof.ergm(), ..., obs.twostage=FALSE){
  if(obs.twostage && control$nsim %% obs.twostage !=0) stop("Number of imputation networks specified by obs.twostage= must divide the nsim control parameter evenly.")
  
  nw <- object$network
  nnets <- length(unique(.peek_vattrv(nw, ".NetworkID")))

  if(is.numeric(subset)) subset <- unwhich(subset, nnets)
  subset <- rep(subset, length.out=nnets)
  remain <- subset

  stats <- stats.obs <- NULL

  message("Constructing simulation model(s).")
  sim_settings <- simulate(object, monitor=NULL, nsim=control$nsim, control=set.control.class("control.simulate.ergm",control), basis=nw, output="stats", response = object$response, ..., do.sim=FALSE)
  if(!is.null(object$constrained.obs)) sim.obs_settings <- simulate(object, monitor=NULL, observational=TRUE, nsim=control$nsim, control=set.control.class("control.simulate.ergm",control), basis=nw, output="stats", response = object$response, ..., do.sim=FALSE)

  prev.remain <- NULL
  cl <- ergm.getCluster(control)
  nthreads <- max(length(cl),1)
  while(any(remain)){
    message("Constructing GOF model.")
    if(any(NVL(prev.remain,FALSE)!=remain))
      pernet.m <- ergm_model(~N(GOF, subset=remain), nw=nw, response = object$response, ...,
                             term.options= modifyList(as.list(object$control$term.options), list(N.compact_stats=FALSE)))
    prev.remain <- remain
    
    nstats <- nparam(pernet.m, canonical=TRUE)/sum(remain)
    
    message("Simulating unconstrained sample.")
    sim <- do.call(simulate, modifyList(sim_settings, list(monitor=pernet.m)))
    sim <- sim[,ncol(sim)-sum(remain)*nstats+seq_len(sum(remain)*nstats),drop=FALSE]

    if(!is.null(object$constrained.obs)){
      sim.obs2 <-
        if(obs.twostage){
          message("Simulating imputed networks.", appendLF=FALSE)
          sim.obs <- list()
          sim.obs.net <- sim_settings$basis
          genseries <- function(){
            for(i in seq_len(ceiling(obs.twostage/nthreads))){
              args <- modifyList(sim_settings,
                                 list(
                                   basis=sim.obs.net,
                                   control=modifyList(sim_settings$control,
                                                      list(parallel=0,MCMC.burnin=if(i==1)sim_settings$control$MCMC.burnin else sim_settings$control$MCMC.interval)),
                                   output="pending_update_network", nsim=1))
              sim.obs.net <- do.call(simulate, args)
              args <- modifyList(sim.obs_settings,
                                 list(basis=sim.obs.net, monitor=pernet.m, nsim=control$nsim/obs.twostage,
                                      control=modifyList(sim.obs_settings$control,
                                                         list(parallel=0))))
              sim.obs[[i]] <- do.call(simulate, args)
              sim.obs[[i]] <- sim.obs[[i]][,ncol(sim.obs[[i]])-sum(remain)*nstats+seq_len(sum(remain)*nstats),drop=FALSE]
              message(".", appendLF=FALSE)
            }
            sim.obs
          }
          sim.obs <- if(!is.null(cl)) unlist(clusterCall(cl, genseries),recursive=FALSE) else genseries()
          message("")
          do.call(rbind, sim.obs)
        }
      # FIXME: Simulations can be rerun only on the networks in the subset.
      message("Simulating constrained sample.")
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
    if(obs.twostage){
      statarray.obs2 <- array(c(sim.obs2), c(control$nsim, nstats, sum(remain)))
      dimnames(statarray.obs2) <- list(Iterations=NULL, Statistic=cn, Network=NULL)
    }else statarray.obs2 <- NULL

    if(is.null(stats)){
      stats <- statarray
      stats.obs <- statarray.obs
      stats.obs2 <- statarray.obs2
    }else{
      stats[,,remain[subset]] <- statarray
      stats.obs[,,remain[subset]] <- statarray.obs
      stats.obs2[,,remain[subset]] <- statarray.obs2
    }

    # Calculate variances for each network and statistic.
    dv <- apply(statarray, 2:3, var) - NVL3(statarray.obs2, apply(., 2:3, .mean_var, obs.twostage), apply(statarray.obs, 2:3, var))
    # If any statistic for the network has negative variance estimate, rerun it.
    remain[remain] <- apply(dv<=0, 2, any)
    if(any(remain)) message(sum(remain), " networks (", paste(which(remain),collapse=", "), ") have bad simulations; rerunning.")
  }
  o <- structure(list(nw=nw, subset=subset, stats=stats, stats.obs=stats.obs, stats.obs2=stats.obs2, obs.twostage=obs.twostage, control=control), class="gofN")
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
  statarray.obs2 <- x$stats.obs2
  control <- x$control
  cn <- dimnames(statarray)[[2]]
  for(i in seq_along(cn)){
    s <- statarray[,i,]
    so <- statarray.obs[,i,]
    so2 <- statarray.obs2[,i,]

    fitted <- colMeans(s)
    svar <- apply(s,2,var)
    sovar <- if(x$obs.twostage) apply(so2,2,.mean_var,x$obs.twostage) else apply(so,2,var)
    resid <- (colMeans(so)-colMeans(s))/sqrt(svar-sovar)

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
