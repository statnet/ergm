#  File R/ergm.stepping.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################

ergm.stepping = function(init, nw, model, initialfit, constraints,
                         control, proposal, proposal.obs, 
                         verbose=FALSE, ...){

  .Deprecate_once(msg="All improvements introduced by the Stepping algorithm have been integrated into the much more flexible MCMLE implementation, so ergm.stepping() is deprecated and will be removed in a future version.")
  control <- remap_algorithm_MCMC_controls(control, "Step")

  #   preliminary, to set up structure. 
  nw.orig <- nw
  asyse=init-init
  mc.se=1+0.05*asyse
  mle.lik=initialfit$mle.lik
  theta.original=init
  
  ## Prepare the output structure:
  obsstats <- summary(model, nw)  # Observed statistics
  init <- init  # beginning parameter value
  samples <- list()  # matrices of sampled network statistics
  sampmeans <- list() # vectors of column means of stats matrices
  xi <- list() # "new obsstats" values, somewhere between obsstats and sampmeans
  eta <- list() # MLE using lognormal approximation and "new obsstats"
  gamma <- list() # factor controlling convex combo: # xi=gamma*obsstats + (1-gamma)*sampmeans  
  
  iter <- 0
  eta[[1]] <- init
  finished <- FALSE
  countdown <- 2
  while (!finished) { # Iterate until gamma==1
    iter=iter+1
    ## Generate an mcmc sample from the probability distribution determined by orig.mle
    samples[[iter]]=simulate(model, nsim=control$MCMC.samplesize, basis=nw,
                                     coef=eta[[iter]], output="stats",
                                     constraints=constraints, 
                                     control=set.control.class("control.simulate.formula",control), ...)
    sampmeans[[iter]]=colMeans(samples[[iter]])
    
    hi <- control$Step.gridsize  # Goal: Let gamma be largest possible multiple of .01
    lo <- gamm <- 0
    message("Iteration #",iter, ". ",appendLF=FALSE)
    if (verbose) {
      message("Current canonical parameter:")
      message_print(eta[[iter]])
    }
    while (hi-lo>1 || hi > gamm) {
      gamm<-ceiling((hi+lo)/2)
      gamma[[iter]] = gamm/control$Step.gridsize
      xi[[iter]] = gamma[[iter]]*obsstats  + (1-gamma[[iter]])*sampmeans[[iter]]
      inCH=is.inCH(1.05*gamma[[iter]]*obsstats  + 
                   (1 - 1.05*gamma[[iter]])*sampmeans[[iter]],
                   samples[[iter]])
      if (inCH) {
        lo <-gamm
        # If 1/gridsize is not small enough, function gives warning message
      } else {
        hi <- gamm
        if (gamm==1) {
          warning("gamma=", 1/control$Step.gridsize, " still not small enough to ",
                  "stay in convex hull.  A larger gridsize than ",
                  control$Step.gridsize, " might help.")
        }
      }
    }
    if (!inCH && gamm>1) {# Last attempt was a fail, so decrease gamm by one
      gamm <- gamm-1
      gamma[[iter]] <- gamm/control$Step.gridsize
      xi[[iter]] <- gamma[[iter]]*obsstats  + (1-gamma[[iter]])*sampmeans[[iter]]
    }
  
    # Now we have found a gamm that moves xi inside the convex hull, 
    # a bit away from the boundary.  This is described 
    # in Hummel, Hunter, Handcock (2011, JCGS).
    if (gamm == control$Step.gridsize) {
      if (verbose) 
        message("Observed stats are well inside the convex hull.")
      countdown <- countdown - 1
    } else {
      countdown <- 2 #  Reset countdown
    }
    # We'd like to have gamma==1 for 2 consecutive iterations before
    # we declare that we're finished.
    finished = (countdown==0) || (iter >= control$Step.maxit)
    
    # When the stepped xi is in the convex hull (but not on the boundary), find the MLE for gyobs=xi
    message("  Trying gamma=", gamma[[iter]],"")
    #' @importFrom utils flush.console
    flush.console()
    
    ############# PLOTS print if VERBOSE=TRUE #################
    if(verbose){
      # Take a look at obsstats (in red dot) and the "new obsstats" (in green triangle):
      # par(mgp = c(2,.8,0), mar = .1+c(3,3,3,1)) ## How do I put more margin at the top?
      #' @importFrom graphics par
      par(ask=TRUE)
      #' @importFrom graphics pairs
      pairs(rbind(samples[[iter]], obsstats, xi[[iter]], sampmeans[[iter]]), 
            col=c(rep(1, nrow(samples[[iter]])), 2, 7, 3), # all black used for JCGS article 
            pch=c(rep(46, nrow(samples[[iter]])), 3, 16, 4),
            cex=c(rep(1, nrow(samples[[iter]])), 1.4, 1, 1.4),
            main=paste("Iteration ", iter, ": ",
                       "gamma = ", gamma[[iter]], sep=""),
            cex.main=1.5, cex.axis=1.5, cex.lab=1.5)
    } #ends if(verbose)
    ############# END PLOTS #################
    # ergm.estimate requires that the simulated stats be "centered" in the sense that the
    # observed statistics would be exactly zero on the same scale.  In this case, the
    # "observed statistics" equal xi[[iter]].
    v<-ergm.estimate(init=eta[[iter]], model=model, 
                     statsmatrices=mcmc.list(as.mcmc(sweep(samples[[iter]], 2, xi[[iter]], '-'))), 
                     nr.maxit=control$MCMLE.NR.maxit,
                     metric=control$MCMLE.metric,
                     verbose=verbose,
                     trace=0,  # suppress 'optim' output
                     estimateonly=TRUE, 
                     #statsmatrix.obs=statsmatrix.obs, 
                     #epsilon=control$epsilon,
                     # nr.reltol=control$MCMLE.NR.reltol,
                     #calc.mcmc.se=control$MCMC.addto.se, hessianflag=control$main.hessian,
                     # method=control$MCMLE.method,
                     ...)
    eta[[iter+1]]<-coef(v)
  }
  message("Now ending with one large sample for MLE. ")
  flush.console()
  iter <- iter+1
  finalsample <- simulate(model, nsim=control$MCMC.samplesize, basis=nw,
                                  coef=eta[[iter]], output="stats", 
                                  constraints=constraints, 
                                  control=set.control.class("control.simulate.formula",control), ...)
  sampmeans[[iter]] <- colMeans(finalsample)
  xi[[iter]] <- obsstats
  v<-ergm.estimate(init=eta[[iter]], model=model, 
                   statsmatrices=mcmc.list(as.mcmc(sweep(finalsample, 2, xi[[iter]], '-'))), 
                   nr.maxit=control$MCMLE.NR.maxit,
                   metric=control$MCMLE.metric,
                   verbose=verbose,
                   trace=0,  # suppress 'optim' output
                   #estimateonly=TRUE,
                   #statsmatrix.obs=statsmatrix.obs, 
                   epsilon=control$epsilon,
                    nr.reltol=control$MCMLE.NR.reltol,
                   calc.mcmc.se=control$MCMC.addto.se, hessianflag=control$main.hessian,
                    method=control$MCMLE.method,
                   ...)
  eta[[iter+1]] <- coef(v)
  
  ############# FINAL PLOT 1 prints if VERBOSE=TRUE #################
  if(verbose){
    pairs(rbind(finalsample, obsstats, sampmeans[[iter]]), 
          col=c(rep(1, nrow(finalsample)), 2, 3),# all black used for JCGS article 
          pch=c(rep(46, nrow(finalsample)),3 ,4 ),
          cex=c(rep(1, nrow(finalsample)), 1.4, 1.4),
          main=paste("Final Stepping Iteration (#", iter, ")", sep=""),# "\n",
          cex.main=1.5, cex.axis=1.5, cex.lab=1.5)
  } #ends if(verbose)
  ############## END PLOT #################    
  #####  final.mle
  mle.lik <- mle.lik + abs(v$loglikelihood)
  v$newnetwork <- nw
  v$burnin <- control$MCMC.burnin
  v$samplesize <- control$MCMC.samplesize
  v$interval <- control$MCMC.interval
  v$network <- nw.orig
  v$newnetwork <- nw
  v$interval <- control$MCMC.interval
  v$theta.original <- theta.original
  v$mplefit <- initialfit
  v$parallel <- control$parallel
  # The following output is sometimes helpful.  It's the 
  # total history of all eta values along with all of the corresponding
  # mean value parameter estimates:
  v$allmeanvals <- t(sapply(sampmeans, function(a)a))
  v$allparamvals <- t(sapply(eta, function(a)a))
  
  if(!v$failure & !any(is.na(coef(v)))){
    asyse <- mc.se
    if(is.null(v$covar)){
      asyse[names(coef(v))] <- suppressWarnings(sqrt(diag(sginv(-v$hessian))))
    }else{
      asyse[names(coef(v))] <- suppressWarnings(sqrt(diag(v$covar)))
    }
  }
  
  v$sample <- ergm.sample.tomcmc(v$sample, control)
  v$etamap <- model$etamap
  v$iterations <- iter

  v
}  # Ends the whole function


#' Transform points represented by rows of `m` such that their
#' centroid is shifted `1-gamma.tr` of the way toward `x` and their
#' spread around their centroid is scaled by a factor of
#' `gamma.scl`. Both `gamma.tr` and `gamma.scl` can be vectors equal
#' in length to the number of columns of `m`.
#' @noRd
.shift_scale_points <- function(m, x, gamma.tr, gamma.scl=gamma.tr){
  m.c <- colMeans(m)
  m <- sweep_cols.matrix(m, m.c, disable_checks=TRUE)
  m.c <- c(m.c*gamma.tr + x*(1-gamma.tr))
  t(t(m) * gamma.scl + m.c)
}

## Given two matrices x1 and x2 with d columns (and any positive
## numbers of rows), find the greatest gamma<=steplength.max s.t., the
## points of x2 shrunk towards the centroid of x1 a factor of gamma,
## are all in the convex hull of x1, as is the centroid of x2 shrunk
## by margin*gamma.
##
## If -1 <= margin < 0, the algorithm "forgives" some amount of either
## centroid or x2 being outside of the convex hull of x1.

## This is a variant of Hummel et al. (2010)'s steplength algorithm
## also usable for missing data MLE.
.Hummel.steplength <- function(x1, x2=NULL, margin=0.05, steplength.max=1, x2.num.max=ceiling(sqrt(ncol(rbind(x1)))), parallel=c("observational","never"), control=NULL, verbose=FALSE){
  parallel <- match.arg(parallel)
  margin <- 1 + margin
  x1 <- rbind(x1); m1 <- rbind(colMeans(x1)); ; n1 <- nrow(x1)
  if(is.function(x2.num.max)) x2.num.max <- x2.num.max(x1)
  if(is.null(x2)){
    m2 <- rbind(rep(0,ncol(x1)))
    parallel <- FALSE
  }else{                                      
    x2 <- rbind(x2)
    m2 <- rbind(colMeans(x2))
    parallel <- parallel == "observational"
  }
  n2 <- nrow(x2)

  # Drop duplicated elements in x1, *and* those elements in x2 that
  # duplicate those in x1.
  d12 <- duplicated(rbind(x1,x2))
  d1 <- d12[seq_len(n1)]
  d2 <- d12[-seq_len(n1)]
  x1 <- x1[!d1,,drop=FALSE]
  x2 <- x2[!d2,,drop=FALSE]
  if(length(x2)==0) x2 <- NULL

  if(verbose>1) message("Eliminating repeated points: ", sum(d1),"/",length(d1), " from target set, ", sum(d2),"/",length(d2)," from test set.")

  cl <- if(parallel && !is.null(control)) ergm.getCluster(control, verbose)

  ## Use PCA to rotate x1 into something numerically stable and drop
  ## unused dimensions, then apply the same affine transformation to
  ## m1 and x2:
  if(nrow(x1)>1){
    ## Center:
    x1m <- colMeans(x1) # note that colMeans(x)!=m1
    x1c <- sweep(x1, 2, x1m, "-")
    ## Rotate x1 onto its principal components, dropping linearly dependent dimensions:
    e <- eigen(crossprod(x1c), symmetric=TRUE)
    Q <- e$vec[,sqrt(pmax(e$val,0)/max(e$val))>sqrt(.Machine$double.eps)*2,drop=FALSE]
    x1cr <- x1c%*%Q # Columns of x1cr are guaranteed to be linearly independent.

    ## Scale x1:
    x1crsd <- pmax(apply(x1cr, 2, sd), sqrt(.Machine$double.eps))
    x1crs <- sweep(x1cr, 2, x1crsd, "/")

    ## Now, apply these operations to m1 and x2:
    m1crs <- sweep(sweep(m1, 2, x1m, "-")%*%Q, 2, x1crsd, "/")
    if(!is.null(x2)) x2crs <- sweep(sweep(x2, 2, x1m, "-")%*%Q, 2, x1crsd, "/")
    m2crs <- sweep(sweep(m2, 2, x1m, "-")%*%Q, 2, x1crsd, "/")
  }else{
    if(is.null(x2)){
      if(isTRUE(all.equal(m1,m2,check.attributes=FALSE))) return(1) else return(0)
    }else{
      if(apply(x2, 1, all.equal, m1, check.attributes=FALSE) %>% map_lgl(isTRUE) %>% all) return(1) else return(0)
    }
  }

  if(!is.null(x2) && nrow(x2crs) > x2.num.max){
    ## If constrained sample size > x2.num.max
    if(verbose>1){message("Using fast and approximate Hummel et al search.")}
    d <- rowSums(sweep(x2crs, 2, m1crs)^2)
    x2crs <- x2crs[order(-d)[1:x2.num.max],,drop=FALSE]
  }

  if(!is.null(cl)){
    # NBs: parRapply() would call shrink_into_CH() for every
    # row. Direct reference to split.data.frame() is necessary here
    # since no matrix method.
    x2crs <- split.data.frame(x2crs, rep_len(seq_along(cl), nrow(x2crs)))
    min(steplength.max, unlist(parallel::parLapply(cl=cl, x2crs, shrink_into_CH, x1crs, verbose=verbose))/margin)
  }else{
    min(steplength.max, shrink_into_CH(if(!is.null(x2)) x2crs else m2crs, x1crs, verbose=verbose)/margin)
  }
}
