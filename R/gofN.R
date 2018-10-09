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
#' @param \dots additional arguments to functions ([simulate.ergm()]
#'   and [summary.ergm_model()] for the constructor, [plot()],
#'   [qqnorm()], and [qqline()] for the plotting method).
#' @param control See [control.gofN.ergm()].
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
#' \item{pearson}{The Pearson residual computed from the above.}
#' 
#' In addition, the following [`attr`]-style attributes are included:
#'
#' \item{nw}{The observed multinetwork object.}
#' 
#' \item{subset}{A logical vector giving the subset of networks that were used.}
#' 
#' \item{control}{Control parameters passed.}
#'
#' @examples
#' data(samplk)
#' monks <- Networks(samplk1, samplk2, samplk3,samplk1, samplk2, samplk3,samplk1, samplk2, samplk3)
#' fit <- ergm(monks~N(~edges))
#' fit.gof <- gofN(fit, GOF=~edges)
#' summary(fit.gof)
#' plot(fit.gof)
#'
#' samplk1[1,]<-NA
#' samplk2[,2]<-NA
#' monks <- Networks(samplk1, samplk2, samplk3,samplk1, samplk2, samplk3,samplk1, samplk2, samplk3)
#' fit <- ergm(monks~N(~edges))
#' fit.gof <- gofN(fit, GOF=~edges)
#' summary(fit.gof)
#' plot(fit.gof)
#' 
#' # Default is good enough in this case, but sometimes, we might want to set it higher. E.g.,
#' \dontrun{
#' fit.gof <- gofN(fit, GOF=~edges, control=control.gofN.ergm(nsim=400))
#' }
#' 
#' @export
gofN <- function(object, GOF, subset=TRUE, control=control.gofN.ergm(), ...){
  check.control.class(c("gofN.ergm"), "gofN")
  if(control$obs.twostage && control$nsim %% control$obs.twostage !=0) stop("Number of imputation networks specified by obs.twostage control parameter must divide the nsim control parameter evenly.")
  
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
  }else control$obs.twostage <- FALSE # Ignore two-stage setting if no observational process.

  prev.remain <- NULL
  cl <- ergm.getCluster(control)
  nthreads <- max(length(cl),1)
  for(attempt in seq_len(control$retry_bad_nets + 1)){
    if(attempt!=1) message(sum(remain), " networks (", paste(which(remain),collapse=", "), ") have bad simulations; rerunning.")
    message("Constructing GOF model.")
    if(any(NVL(prev.remain,FALSE)!=remain))
      pernet.m <- ergm_model(~N(GOF, subset=remain), nw=nw, response = object$response, ...,
                             term.options= .update.list(as.list(object$control$term.options), list(N.compact_stats=FALSE)))
    prev.remain <- remain
    
    nstats <- nparam(pernet.m, canonical=TRUE)/sum(remain)

    # TODO: Simulations can be rerun only on the networks in the subset.
    
    # The two-stage sample, taken marginally, *is* an unconstrained
    # sample.
    if(!control$obs.twostage){
      message("Simulating unconstrained sample.")
      sim <- do.call(simulate, .update.list(sim_settings, list(monitor=pernet.m)))
      sim <- sim[,ncol(sim)-sum(remain)*nstats+seq_len(sum(remain)*nstats),drop=FALSE]
    }

    # TODO: Make this adaptive: start with a small simulation,
    # increase on fail; or perhaps use a pilot sample.
    if(!is.null(object$constrained.obs)){
      sim <-
        if(control$obs.twostage){
          message("Simulating imputed networks.", appendLF=FALSE)
          sim.net <- sim_settings$basis
          genseries <- function(){
            sim <- list()
            for(i in seq_len(ceiling(control$obs.twostage/nthreads))){
              args <- .update.list(sim_settings,
                                 list(
                                   basis=sim.net,
                                   control=.update.list(sim_settings$control,
                                                      list(parallel=0,MCMC.burnin=if(i==1)sim_settings$control$MCMC.burnin else sim_settings$control$MCMC.interval)),
                                   output="pending_update_network", nsim=1))
              sim.net <- do.call(simulate, args)
              args <- .update.list(sim.obs_settings,
                                 list(basis=sim.net, monitor=pernet.m, nsim=control$nsim/control$obs.twostage,
                                      control=.update.list(sim.obs_settings$control,
                                                         list(parallel=0))))
              sim[[i]] <- do.call(simulate, args)
              sim[[i]] <- sim[[i]][,ncol(sim[[i]])-sum(remain)*nstats+seq_len(sum(remain)*nstats),drop=FALSE]
              message(".", appendLF=FALSE)
            }
            sim
          }
          #' @importFrom parallel clusterCall
          sim <- if(!is.null(cl)) unlist(clusterCall(cl, genseries),recursive=FALSE) else genseries()
          message("")
          do.call(rbind, sim)
        }
      message("Simulating constrained sample.")
      sim.obs <- do.call(simulate, .update.list(sim.obs_settings,
                                                list(nsim=if(control$obs.twostage) control$obs.twostage else control$nsim,
                                                     monitor=pernet.m)))
      sim.obs <- sim.obs[,ncol(sim.obs)-sum(remain)*nstats+seq_len(sum(remain)*nstats),drop=FALSE]
    }else{
      sim.obs <- matrix(summary(pernet.m, object$network, response = object$response, ...), nrow(sim), ncol(sim), byrow=TRUE)
    }

    message("Collating the simulations.")
    cn <- colnames(sim)[seq_len(nstats)] %>% sub(".*?:","", .)
    
    statarray <- array(c(sim), c(control$nsim, nstats, sum(remain)))
    dimnames(statarray) <- list(Iterations=NULL, Statistic=cn, Network=NULL)
    statarray.obs <- array(c(sim.obs), c(if(control$obs.twostage) control$obs.twostage else control$nsim, nstats, sum(remain)))
    dimnames(statarray.obs) <- list(Iterations=NULL, Statistic=cn, Network=NULL)

    if(is.null(stats)){
      stats <- statarray
      stats.obs <- statarray.obs
    }else{
      stats[,,remain[subset]] <- statarray
      stats.obs[,,remain[subset]] <- statarray.obs
    }

    # Calculate variances for each network and statistic.
    dv <- apply(statarray, 2:3, var) - if(control$obs.twostage) apply(statarray, 2:3, .mean_var, control$obs.twostage) else apply(statarray.obs, 2:3, var)
    # If any statistic for the network has negative variance estimate, rerun it.
    remain[remain] <- apply(dv<=0, 2, any)
    if(!any(remain)) break;
  }
  if(any(remain))
    stop(sum(remain), " networks (", paste(which(remain),collapse=", "), ") have bad simulations after permitted number of retries. Rerun with higher nsim= control parameter.")

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
        var.obs = if(control$obs.twostage) apply(s, 2, .mean_var, control$obs.twostage) else apply(so, 2, var),
        pearson = (observed-fitted)/sqrt(var-var.obs)
      )
  }), cn)

  structure(o, nw=nw, subset=subset, control=control, class="gofN")
}

#' @describeIn gofN A plotting method, making residual and scale-location plots.
#'
#' @param x a [`gofN`] object.
#' @param against vector of values, network attribute, or a formula whose RHS gives an expression in terms of network attributes to plot against; if `NULL` (default), plots against fitted values.
#' @param col,pch,cex vector of values (wrapped in [I()]), network attribute, or a formula whose RHS gives an expression in terms of network attributes to plot against.
#' @param which which to plot (`1` for residuals plot, `2` for \eqn{\sqrt{|R_i|}}{sqrt(|R_i|)} scale plot, and `3` for normal quantile-quantile plot).
#' @param ask whether the user should be prompted between the plots.
#' 
#' @importFrom grDevices dev.interactive devAskNewPage
#' @importFrom graphics abline panel.smooth plot
#' @importFrom methods is
#' @export
plot.gofN <- function(x, against=NULL, which=1:2, col=1, pch=1, cex=1, ..., ask = length(which)>1 && dev.interactive(TRUE)){
  if(ask){
    prev.ask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(prev.ask))
  }

  if(any(sapply(list(against, col, pch, cex),
                function(x) is.character(x) || is(x,"formula"))))
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
      qqline(summ$pearson, ...)
    }
    
  }
}

#' @describeIn gofN A simple summary function.
#' @param by a numeric or character vector, or a formula whose RHS gives an expression in terms of network attributes, used as a grouping variable for summarizing the values.
#' @export
summary.gofN <- function(object, by=NULL, ...){
  cns <- names(object)
  if(is.null(by)){
    list(`Observed/Imputed values` = object %>% map("observed") %>% as_tibble %>% summary,
         `Fitted values` = object %>% map("fitted") %>% as_tibble %>% summary,
         `Pearson residuals`  = object %>% map("pearson") %>% as_tibble %>% summary,
         `Variance of Pearson residuals` = object %>% map("pearson") %>% map(var),
         `Std. dev. of Pearson residuals` = object %>% map("pearson") %>% map(sd))
  }else{
    if(is(by,"formula"))
      nattrs <- get_multinet_nattr_tibble(attr(object,"nw"))[attr(object,"subset"),]
    
    byname <- switch(class(by),
                     formula = despace(deparse(by[[length(by)]])),
                     despace(deparse(substitute(by))))
    byval <- switch(class(by),
                    formula = eval(by[[length(by)]], envir = nattrs, enclos = environment(by)),
                    by)

    list(`Observed/Imputed values` = object %>% map("observed") %>% as_tibble %>% split(byval) %>% map(summary),
         `Fitted values` = object %>% map("fitted") %>% as_tibble %>% split(byval) %>% map(summary),
         `Pearson residuals`  = object %>% map("pearson") %>% as_tibble %>% split(byval) %>% map(summary),
         `Variance of Pearson residuals` = object %>% map("pearson") %>% as_tibble %>% split(byval) %>% map(var),
         `Std. dev. of Pearson residuals` = object %>% map("pearson") %>% as_tibble %>% split(byval) %>% map(~apply(.,2,sd)))
  }
}

#' Auxiliary for Controlling Multinetwork ERGM Linear Goodness-of-Fit Evaluation
#'
#' @inherit control.gof.ergm description
#' @inheritParams control.gof.ergm
#' @param obs.twostage Either `FALSE` or an integer. This parameter
#'   only has an effect if the network has missing data or
#'   observational process. For such networks, evaluating the Pearson
#'   residual requires simulating the expected value of the
#'   conditional variance under the observation process. If `FALSE`,
#'   the simulation is performed conditional on the observed
#'   network. However, a more accurate estimate can be obtained via a
#'   two-stage process: \enumerate{
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
#' to simulate from, which should divide the [control.gofN.ergm()]'s
#' `nsim` argument evenly.
#' 
#' @param retry_bad_nets This parameter only has an effect if the
#'   network has missing data or observational process. It gives the
#'   number of times to retry simulating networks for which the
#'   estimated constrained variance is higher than unconstrained. Note
#'   that setting it `>0` is likely to bias estimates: the simulation
#'   should instead be rerun with a larger `nsim`.
#' 
#' @description `control.gofN.ergm` (or its alias, `control.gofN`) is
#'   intended to be used with [gofN()] specifically and will "inherit"
#'   as many control parameters from [`ergm`] fit as possible().
#'  
#' @export control.gofN.ergm
control.gofN.ergm<-function(nsim=100,
                            obs.twostage=nsim/2,
                            retry_bad_nets=0,
                       MCMC.burnin=NULL,
                       MCMC.interval=NULL,
                       MCMC.prop.weights=NULL,
                       MCMC.prop.args=NULL,
                       
                       MCMC.init.maxedges=NULL,
                       MCMC.packagenames=NULL,
                       
                       MCMC.runtime.traceplot=FALSE,
                       network.output="network",
                       
                       seed=NULL,
                       parallel=0,
                       parallel.type=NULL,
                       parallel.version.check=TRUE,
                       parallel.inherit.MT=FALSE){
  control<-list()
  for(arg in names(formals(sys.function())))
    control[arg]<-list(get(arg))

  set.control.class("control.gofN.ergm")
}

#' @rdname control.gofN.ergm
#' @export control.gofN
control.gofN <- control.gofN.ergm
