ergm.mainfitloop <- function(theta0, nw, model, Clist,
                             initialfit, 
                             MCMCparams, 
                             MHproposal, MHproposal.miss,
                             verbose=FALSE,
                             epsilon=1e-10,
                             estimate=TRUE, ...) {
  iteration <- 1
  nw.orig <- network.copy(nw)
  asyse=theta0-theta0
  mc.se=1+0.05*asyse
  mle.lik=initialfit$mle.lik
  v=list(coef=theta0)
  theta.original=theta0
  thetaprior=theta0-theta0

  stats <- matrix(0,ncol=Clist$nstats,nrow=MCMCparams$samplesize)
  stats[1,] <- Clist$obs - Clist$meanstats
  MCMCparams$stats <- stats
  MCMCparams$meanstats <- Clist$meanstats
  if(network.naedgecount(nw) > 0){
    MCMCparams.miss <- MCMCparams
    if(!is.null(MCMCparams$miss.MCMCsamplesize)){
      MCMCparams.miss$MCMCsamplesize <- MCMCparams$miss.MCMCsamplesize
    }
    if(!is.null(MCMCparams$miss.interval)){
      MCMCparams.miss$interval <- MCMCparams$miss.interval
    }
    if(!is.null(MCMCparams$miss.burnin)){
      MCMCparams.miss$burnin <- MCMCparams$miss.burnin
    }
  }
  while(iteration <= MCMCparams$maxit){
    thetaprior <- theta0
    theta0 <- v$coef
    eta0 <- ergm.eta(theta0, model$etamap)
    if(verbose){
      cat("Iteration ",iteration," of at most ", MCMCparams$maxit,
          " with parameter: \n", sep="")
      print(theta0)
    }else{
      cat("Iteration ",iteration," of at most ", MCMCparams$maxit,": \n",sep="")
    }
    z <- ergm.getMCMCsample.parallel(nw, model, MHproposal, eta0, MCMCparams, verbose)
    statsmatrix <- z$statsmatrix
    v$sample <- statsmatrix
    nw.returned <- network.copy(z$newnetwork)
    if(network.naedgecount(nw) > 0){
      z.miss <- ergm.getMCMCsample.parallel(nw, model, MHproposal.miss, eta0, MCMCparams.miss, verbose)
      statsmatrix.miss <-z.miss$statsmatrix
      nw.miss.returned <- network.copy(z.miss$newnetwork)
      if(verbose){cat("Back from constrained MCMC...\n")}
    }else{
      statsmatrix.miss <- NULL
      if(verbose){cat("Back from unconstrained MCMC...\n")}
    }
    if(MCMCparams$sequential & network.naedgecount(nw) == 0){
      nw <- nw.returned
      nw.obs <- summary(model$formula, basis=nw)
      namesmatch <- match(names(MCMCparams$meanstats), names(nw.obs))
      MCMCparams$stats[1,] <- nw.obs[namesmatch]-MCMCparams$meanstats
    }
#
    iteration <- iteration + 1
    if(z$nedges >= 50000-1 || ergm.checkdegeneracy(statsmatrix, statsmatrix.miss, verbose=verbose)){
     if(iteration <= MCMCparams$maxit){
      cat(paste("The MCMC sampler is producing degenerate samples.\n",
                "Try starting the algorithm at an alternative model\n",
                "(That is, changing the 'theta0' argument).\n",
                "I am trying something simple...\n",
                "The current theta0 is:\n"))
                print(theta0)
      v$coef <- 0.9*theta0
      next
     }else{
      cat(paste("The MCMC sampler is producing degenerate samples.\n",
                "Try starting the algorithm at an alternative model\n",
                "(That is, changing the 'theta0' argument).\n",
                "The current theta0 is:\n"))
              print(theta0)
      v$coef <- theta0
      return(structure (v, class="ergm"))
     }
    }
    if(verbose){
      cat(paste("The density of the returned network is",
                network.density(nw.returned),"\n"))
      cat(paste("The density of the original network is",
                network.density(nw.orig),"\n"))
      cat("Summary of simulation, relative to observed network:\n")
      b = apply(statsmatrix,2,summary.statsmatrix.ergm)
      print(b,scipen=6)
      degreedist(nw.returned)
      cat("Meanstats of simulation, relative to observed network:\n")
      print(summary(model$formula, basis=nw.returned)-Clist$meanstats)
      if(network.naedgecount(nw) > 0){
       cat("Summary of simulation, relative to missing network:\n")
        a = apply(statsmatrix.miss,2,summary.statsmatrix.ergm)[4,]
        b = sweep(apply(statsmatrix,2,summary.statsmatrix.ergm),2,a,"-")
        print(b,scipen=6)
       degreedist(nw.miss.returned)
       cat("Meanstats of simulation, relative to missing network:\n")
       print(summary(model$formula, basis=nw.miss.returned)-Clist$meanstats)
       nw.returned <- network.copy(nw.miss.returned)
      }
    }
    if(verbose){cat("Calling optimization routines...\n")}

###### old statseval begins here.
    
  if(verbose) {cat("statseval:  epsilon value is ",MCMCparams$epsilon,"\n")}
  if(!estimate){
    if(verbose){cat("Skipping optimization routines...\n")}
    l <- list(coef=theta0, mc.se=rep(NA,length=length(theta0)),
              sample=statsmatrix, sample.miss=statsmatrix.miss,
              iterations=1, MCMCtheta=theta0,
              loglikelihood=NA, #mcmcloglik=NULL, 
              mle.lik=NULL,
              gradient=rep(NA,length=length(theta0)), #acf=NULL,
              samplesize=MCMCparams$samplesize, failure=TRUE,
              newnetwork = nw.returned)
    return(structure (l, class="ergm"))
  }
  statsmatrix.0 <- statsmatrix
  statsmatrix.0.miss <- statsmatrix.miss
  if(MCMCparams$steplength=="adaptive"){
   if(verbose){cat("Calling adaptive MCMLE Optimization...\n")}
   adaptive.steplength <- 2
   statsmean <- apply(statsmatrix.0,2,mean)
   v <- list(loglikelihood=MCMCparams$adaptive.trustregion*2)
   while(v$loglikelihood > MCMCparams$adaptive.trustregion){
    adaptive.steplength <- adaptive.steplength / 2
    if(!is.null(statsmatrix.0.miss)){
      statsmatrix.miss <- statsmatrix.0.miss*adaptive.steplength+statsmatrix.0*(1-adaptive.steplength)
    }else{
      statsmatrix <- sweep(statsmatrix.0,2,(1-adaptive.steplength)*statsmean,"-")
    }
    if(verbose){cat(paste("Using Newton-Raphson Step with step length",adaptive.steplength,"...\n"))}
#
#   If not the last iteration do not compute all the extraneous
#   statistics that are not needed until output
#
     v<-ergm.estimate(theta0=theta0, model=model,
                      statsmatrix=statsmatrix, 
                      statsmatrix.miss=statsmatrix.miss, 
                      epsilon=MCMCparams$epsilon,
                      nr.maxit=MCMCparams$nr.maxit,
                      nr.reltol=MCMCparams$nr.reltol,
                      calc.mcmc.se=MCMCparams$calc.mcmc.se, hessianflag=MCMCparams$hessian,
                      trustregion=MCMCparams$trustregion, method=MCMCparams$method,
                      metric=MCMCparams$metric,
                      compress=MCMCparams$compress, verbose=verbose,
                      estimateonly=TRUE)
   }
   if(v$loglikelihood < MCMCparams$trustregion-0.001){
     current.scipen <- options()$scipen
    options(scipen=3)
    cat("the log-likelihood improved by",
        format.pval(v$loglikelihood,digits=4,eps=1e-4),"\n")
    options(scipen=current.scipen)
   }else{
    cat("the log-likelihood did not improve.\n")
   }
   if((adaptive.steplength==1) && (v$loglikelihood < MCMCparams$adaptive.epsilon) ){break}
  }else{
    if(verbose){cat("Calling MCMLE Optimization...\n")}
    statsmean <- apply(statsmatrix.0,2,mean)
    if(!is.null(statsmatrix.0.miss)){
      statsmatrix.miss <- statsmatrix.0.miss*MCMCparams$steplength+statsmatrix.0*(1-MCMCparams$steplength)
    }else{
      statsmatrix <- sweep(statsmatrix.0,2,(1-MCMCparams$steplength)*statsmean,"-")
    }
    if(verbose){cat(paste("Using Newton-Raphson Step with step length ",MCMCparams$steplength," ...\n"))}
    v<-ergm.estimate(theta0=theta0, model=model,
                     statsmatrix=statsmatrix, 
                     statsmatrix.miss=statsmatrix.miss, 
                     epsilon=MCMCparams$epsilon,
                     nr.maxit=MCMCparams$nr.maxit,
                     nr.reltol=MCMCparams$nr.reltol,
                     calc.mcmc.se=MCMCparams$calc.mcmc.se, hessianflag=MCMCparams$hessian,
                     trustregion=MCMCparams$trustregion, method=MCMCparams$method,
                     metric=MCMCparams$metric,
                     compress=MCMCparams$compress, verbose=verbose,
                     estimateonly=TRUE)
   if(v$loglikelihood < MCMCparams$trustregion-0.001){
     current.scipen <- options()$scipen
    options(scipen=3)
    cat("the log-likelihood improved by",
        format.pval(v$loglikelihood,digits=4,eps=1e-4),"\n")
    options(scipen=current.scipen)
   }else{
    cat("the log-likelihood did not improve.\n")
   }
   if((MCMCparams$steplength==1) && (v$loglikelihood < MCMCparams$adaptive.epsilon) ){break}
  }
#
# End main loop
#
  }
#
# This is the last iteration, so compute all the extraneous
# statistics that are not needed until output
#
  v<-ergm.estimate(theta0=theta0, model=model,
                   statsmatrix=statsmatrix, statsmatrix.miss=statsmatrix.miss, 
                   epsilon=MCMCparams$epsilon,
                   nr.maxit=MCMCparams$nr.maxit,
                   nr.reltol=MCMCparams$nr.reltol,
                   calc.mcmc.se=MCMCparams$calc.mcmc.se, 
		   hessianflag=MCMCparams$hessian,
                   trustregion=MCMCparams$trustregion, method=MCMCparams$method,
                   metric=MCMCparams$metric,
                   compress=MCMCparams$compress, verbose=verbose)
#
#    Reform the output
#
# check for degeneracy?  See ergm.checkdegeneracy
# l$bounddeg <- list(condAllDegExact=condAllDegExact, maxout=maxout,
#                    maxin=maxin, minout=minout,
#                    minin=minin, attribs=attribs)
#   This is the MCMC estimate as in the Hunter&Handcock paper
#  l$mcmcloglik <- l$mcmcloglik - Clist$ndyads*log(2)
#   This is the MPLE value estimate plus the ratio
#   set in the last line of "ergm.R" 
#   v$mle.lik <- -glm.fit$deviance/2 + v$loglikelihood
  mle.lik <- mle.lik + abs(v$loglikelihood)

# v$newnetwork <- network.copy(z$newnetwork)
# v$newnetwork <- nw

###### old ergm.statseval ends here    
    
    v$burnin <- MCMCparams$burnin
    v$samplesize <- MCMCparams$samplesize
    v$interval <- MCMCparams$interval
    v$network <- network.copy(nw.orig)
    v$newnetwork <- network.copy(nw.returned)
    v$interval <- MCMCparams$interval
    v$theta.original <- theta.original
    v$mplefit <- initialfit
    v$parallel <- MCMCparams$parallel
         
    if(!v$failure & !any(is.na(v$coef))){
#     asyse <- sqrt(diag(robust.inverse(-v$hessian)))
#     asyse <- try(sqrt(diag(robust.inverse(cov(statsmatrix)))))
      asyse <- mc.se
      options(warn=-1)
#     options(warn=2)
      if(is.null(v$covar)){
        asyse[names(v$coef)] <- sqrt(diag(robust.inverse(-v$hessian)))
      }else{
        asyse[names(v$coef)] <- sqrt(diag(v$covar))
      }
      options(warn=0)
    }

  endrun <- MCMCparams$burnin+MCMCparams$interval*(MCMCparams$samplesize-1)
  attr(v$sample, "mcpar") <- c(MCMCparams$burnin+1, endrun, MCMCparams$interval)
  attr(v$sample, "class") <- "mcmc"
# v$null.deviance <- v$mplefit$glm$null.deviance
  v$null.deviance <- 2*network.dyadcount(nw.orig)*log(2)
# v$mle.lik <- -v$mplefit$glm$deviance/2 + v$loglikelihood
# v$mle.lik <- -v$null.deviance/2 + v$loglikelihood
# v$mle.lik <- -initialfit$mle.lik/2 + v$loglikelihood
  v$mle.lik <- mle.lik
# This next is the right one
# v$mcmcloglik <- v$mcmcloglik - network.dyadcount(nw.orig)*log(2)
  v$etamap <- model$etamap
# if(!is.na(v$mle.lik) && v$mle.lik>0){
#   v$mle.lik <- -v$mplefit$glm$deviance/2 
# }  #I'm really not sure about this mle.lik value.
  v
}
