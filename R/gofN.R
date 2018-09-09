.mean_var <- function(x, ng){
  split(x, (seq_along(x)-1L)%/%(length(x)/ng)) %>%
    map_dbl(var) %>%
    mean()
}

.update.list <- function(l, v){
  l[names(v)]<-v
  l
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
#' @return An object of class `gofN`: a named list containing a list
#'   for every statistic in the specified `GOF` formula with the
#'   following elements vectors of length equal to the number of
#'   subnetworks:
#'
#' \item{observed}{For completely observed networks, their value of
#' the statistic. For partially observed networks, the expected value
#' of their imputations under the model.}
#'
#' \item{fitted}{Expected value of the statistic under the model.}
#'
#' \item{var}{Variance of the statistic under the model.}
#'
#' \item{var.obs}{Conditional variance under imputation statistic.}
#' 
#' In addition, the following [`attr`]-style attributes are included:
#'
#' \item{nw}{The observed multinetwork object.}
#' 
#' \item{subset}{A logical vector giving the subset of networks that were used.}
#' 
#' \item{obs.twostage}{Number of runs in the first stage, or `FALSE`.}
#'
#' \item{control}{Control parameters passed.}
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
  
  if(!is.null(object$constrained.obs)){
    sim.obs_settings <- simulate(object, monitor=NULL, observational=TRUE, nsim=control$nsim, control=set.control.class("control.simulate.ergm",control), basis=nw, output="stats", response = object$response, ..., do.sim=FALSE)
  }else obs.twostage <- FALSE # Ignore two-stage setting if no observational process.

  prev.remain <- NULL
  cl <- ergm.getCluster(control)
  nthreads <- max(length(cl),1)
  while(any(remain)){
    message("Constructing GOF model.")
    if(any(NVL(prev.remain,FALSE)!=remain))
      pernet.m <- ergm_model(~N(GOF, subset=remain), nw=nw, response = object$response, ...,
                             term.options= .update.list(as.list(object$control$term.options), list(N.compact_stats=FALSE)))
    prev.remain <- remain
    
    nstats <- nparam(pernet.m, canonical=TRUE)/sum(remain)

    # FIXME: Simulations can be rerun only on the networks in the subset.
    
    # The two-stage sample, taken marginally, *is* an unconstrained
    # sample.
    if(!obs.twostage){
      message("Simulating unconstrained sample.")
      sim <- do.call(simulate, .update.list(sim_settings, list(monitor=pernet.m)))
      sim <- sim[,ncol(sim)-sum(remain)*nstats+seq_len(sum(remain)*nstats),drop=FALSE]
    }

    if(!is.null(object$constrained.obs)){
      sim <-
        if(obs.twostage){
          message("Simulating imputed networks.", appendLF=FALSE)
          sim.net <- sim_settings$basis
          genseries <- function(){
            sim <- list()
            for(i in seq_len(ceiling(obs.twostage/nthreads))){
              args <- .update.list(sim_settings,
                                 list(
                                   basis=sim.net,
                                   control=.update.list(sim_settings$control,
                                                      list(parallel=0,MCMC.burnin=if(i==1)sim_settings$control$MCMC.burnin else sim_settings$control$MCMC.interval)),
                                   output="pending_update_network", nsim=1))
              sim.net <- do.call(simulate, args)
              args <- .update.list(sim.obs_settings,
                                 list(basis=sim.net, monitor=pernet.m, nsim=control$nsim/obs.twostage,
                                      control=.update.list(sim.obs_settings$control,
                                                         list(parallel=0))))
              sim[[i]] <- do.call(simulate, args)
              sim[[i]] <- sim[[i]][,ncol(sim[[i]])-sum(remain)*nstats+seq_len(sum(remain)*nstats),drop=FALSE]
              message(".", appendLF=FALSE)
            }
            sim
          }
          sim <- if(!is.null(cl)) unlist(clusterCall(cl, genseries),recursive=FALSE) else genseries()
          message("")
          do.call(rbind, sim)
        }
      # NOTE: It might be worth using a smaller sample size here, perhaps obs.twostage.
      message("Simulating constrained sample.")
      sim.obs <- do.call(simulate, .update.list(sim.obs_settings, list(monitor=pernet.m)))
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
    dv <- apply(statarray, 2:3, var) - if(obs.twostage) apply(statarray, 2:3, .mean_var, obs.twostage) else apply(statarray.obs, 2:3, var)
    # If any statistic for the network has negative variance estimate, rerun it.
    remain[remain] <- apply(dv<=0, 2, any)
    if(any(remain)) message(sum(remain), " networks (", paste(which(remain),collapse=", "), ") have bad simulations; rerunning.")
  }

  message("Summarising.")
  o <- setNames(lapply(seq_along(cn), function(i){
    s <- stats[,i,]
    so <- stats.obs[,i,]
    #' @importFrom tibble lst
    l <-
      lst(
        observed = colMeans(so),
        fitted = colMeans(s),
        var = apply(s,2,var),
        var.obs = if(obs.twostage) apply(s, 2, .mean_var, obs.twostage) else apply(so, 2, var),
        pearson = (observed-fitted)/sqrt(var-var.obs)
      )
  }), cn)

  structure(o, nw=nw, subset=subset, obs.twostage=obs.twostage, control=control, class="gofN")
}

#' @describeIn gofN A plotting method, making residual and scale-location plots.
#' 
#' @param against vector of values, network attribute, or a formula whose RHS gives an expression in terms of network attributes to plot against; if `NULL` (default), plots against fitted values.
#' @param col,pch,cex vector of values (wrapped in [I()]), network attribute, or a formula whose RHS gives an expression in terms of network attributes to plot against.
#' @param which which to plot (`1` for residuals plot, `2` for \eqn{\sqrt{|R_i|}}{sqrt(|R_i|)} scale plot, and `3` for normal quantile-quantile plot).
#' 
#' @export
plot.gofN <- function(x, against=NULL, which=1:2, col=1, pch=1, cex=1, ..., ask = dev.interactive()){
  if(ask){
    prev.ask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(prev.ask))
  }

  nattrs <- get_multinet_nattr_tibble(attr(x,"nw"))[attr(x,"subset"),]
  
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
  
  for(name in names(x)){
    summ <- x[[name]]
    
    if(1L %in% which){
      plot(NVL(againstval,summ$fitted), summ$pearson, col=col, pch=pch, cex=cex,..., main = paste("Residuals vs. Fitted for", sQuote(name)), xlab=againstname, ylab="Pearson residual",type="n")
      panel.smooth(NVL(againstval,summ$fitted), summ$pearson, col=col, pch=pch, cex=cex, ...)
      abline(h=0, lty=3, col="gray")
    }
    
    if(2L %in% which){
      plot(NVL(againstval,summ$fitted), sqrt(abs(summ$pearson)), col=col, pch=pch, cex=cex,..., main = paste("Scale-location plot for", sQuote(name)), xlab=againstname, ylab=expression(sqrt("|Pearson residual|")), type="n")
      panel.smooth(NVL(againstval,summ$fitted), sqrt(abs(summ$pearson)), col=col, pch=pch, cex=cex, ...)
      abline(h=0, lty=3, col="gray")
    }

    if(3L %in% which){
      qqnorm(summ$pearson, col=col, pch=pch, cex=cex,..., main = paste("Normal Q-Q for", sQuote(name)), xlab=againstname, ylab=expression(sqrt("|Pearson residual|")))
      qqline(summ$pearson)
    }
    
  }
}
