ergm.mainfitloop <- function(theta0, nw, model, Clist,
                             initialfit, 
                             MCMCparams, 
                             MHproposal, 
                             verbose=FALSE,
                             epsilon=1e-10,
                             sequential=TRUE,
                             estimate=TRUE, ...) {
  iteration <- 1
  nw.orig <- nw
#
  asyse=theta0-theta0
  mc.se=1+0.05*asyse
  mle.lik=initialfit$mle.lik
  v=list(coef=theta0)
  theta.original=theta0
  thetaprior=theta0-theta0

  stats <- matrix(0,ncol=Clist$nparam,nrow=MCMCparams$samplesize)
  stats[1,] <- Clist$obs - Clist$meanstats
# stats[,]<-  rep(Clist$obs - Clist$meanstats,rep(nrow(stats),Clist$nparam))
  MCMCparams$stats <- stats
  MCMCparams$meanstats <- Clist$meanstats

#  while(any(mcmc.precision*asyse < mc.se, na.rm=TRUE) && iteration <= maxit){
  while(iteration <= MCMCparams$maxit){
    thetaprior <- theta0
    theta0 <- v$coef
    eta0 <- ergm.eta(theta0, model$etamap)
    if(verbose){
     cat("Iteration ",iteration," of at most ", MCMCparams$maxit,
         " with parameter: \n", sep="")
      print(theta0)
    }else{
     cat("Iteration ",iteration," of at most ", MCMCparams$maxit,": ",sep="")
    }
    z <- ergm.getMCMCsample(nw, model, MHproposal, eta0, MCMCparams, verbose)
    statsmatrix=z$statsmatrix
    statsmatrix.miss <- NULL
    if(verbose){cat("Back from unconstrained MCMC...\n")}
    if(sequential & MCMCparams$Clist.miss$nedges == 0){
      nw <- z$newnetwork
      nw.obs <- summary(model$formula, basis=nw)
      namesmatch <- match(names(MCMCparams$meanstats), names(nw.obs))
      MCMCparams$stats[1,] <- nw.obs[namesmatch]-MCMCparams$meanstats
    }
#
    iteration <- iteration + 1
    if(MCMCparams$steplength<1 && iteration < MCMCparams$maxit ){
        statsmean <- apply(statsmatrix,2,mean)
        statsmatrix <- sweep(statsmatrix,2,(1-MCMCparams$steplength)*statsmean,"-")
    }
    if(z$nedges >= 50000-1 || ergm.checkdegeneracy(statsmatrix, statsmatrix.miss, verbose=verbose)){
     if(iteration <= MCMCparams$maxit){
      cat(paste("The MCMC sampler is producing degenerate samples.\n",
                "Try starting the algorithm at an alternative model\n",
                "(That is, changing the 'theta0' argument).\n",
                "I am trying something simple...\n",
                "The current theta0 is:\n"))
                print(theta0)
#     shrink <- ergm(nw ~ edges)$coef
#     theta0 <- 0.4*theta0
#     theta0["edges"] <- theta0["edges"] + 0.6*shrink
#     v <- list(coef=theta0)
#     v$coef <- theta0
      v$coef <- 0.9*theta0
      next
     }else{
      cat(paste("The MCMC sampler is producing degenerate samples.\n",
                "Try starting the algorithm at an alternative model\n",
                "(That is, changing the 'theta0' argument).\n",
                "The current theta0 is:\n"))
              print(theta0)
#     v <- list(coef=theta0)
      v$coef <- theta0
      return(structure (v, class="ergm"))
     }
    }
    if(verbose){
      cat(paste("The density of the returned network is",
                network.density(nw),"\n"))
      cat(paste("The density of the original network is",
                network.density(nw.orig),"\n"))
      cat("Summary of simulation, relative to observed network:\n")
      print(apply(statsmatrix,2,summary.statsmatrix.ergm),scipen=6)
      degreedist(nw)
    }
    if(verbose){cat("Calling optimization routines...\n")}

###### old statseval begins here.
    
  if(verbose) {cat("statseval:  epsilon value is ",MCMCparams$epsilon,"\n")}
  if(!estimate){
    if(verbose){cat("Skipping optimization routines...\n")}
    l <- list(coef=theta0, mc.se=rep(NA,length=length(theta0)),
              sample=statsmatrix, # sample.miss=statsmatrix.miss,
              iterations=1, MCMCtheta=theta0,
              loglikelihood=NA, #mcmcloglik=NULL, 
              mle.lik=NULL,
              gradient=rep(NA,length=length(theta0)), #acf=NULL,
              samplesize=MCMCparams$samplesize, failure=TRUE,
              newnetwork = nw)
    return(structure (l, class="ergm"))
  }
  if(verbose){cat("Calling MCMLE Optimization...\n")}
  if(verbose){cat("Using Newton-Raphson Step ...\n")}
#
# If not the last iteration do not compute all the extraneous
# statistics that are not needed until output
#
  if(iteration <= MCMCparams$maxit){
   v<-ergm.estimate(theta0=theta0, model=model,
                    statsmatrix=statsmatrix, 
                    statsmatrix.miss=statsmatrix.miss, 
                    epsilon=MCMCparams$epsilon,
                    nr.maxit=MCMCparams$nr.maxit, calc.mcmc.se=MCMCparams$calc.mcmc.se, hessian=MCMCparams$hessian,
                    trustregion=MCMCparams$trustregion, method=MCMCparams$method, metric="Likelihood",
                    compress=MCMCparams$compress, verbose=verbose, estimateonly=TRUE)
  }
#
# End main loop
#
  }
# This is the last iteration, so compute all the extraneous
# statistics that are not needed until output
#
  v<-ergm.estimate(theta0=theta0, model=model,
                   statsmatrix=statsmatrix, statsmatrix.miss=statsmatrix.miss, 
                   epsilon=MCMCparams$epsilon,
                   nr.maxit=MCMCparams$nr.maxit, calc.mcmc.se=MCMCparams$calc.mcmc.se, hessian=MCMCparams$hessian,
                   trustregion=MCMCparams$trustregion, method=MCMCparams$method, metric="Likelihood",
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

  v$newnetwork <- nw

###### old ergm.statseval ends here    
    
    v$burnin <- MCMCparams$burnin
    v$samplesize <- MCMCparams$samplesize
    v$interval <- MCMCparams$interval
    v$network <- nw.orig
    v$newnetwork <- nw
    v$interval <- MCMCparams$interval
    v$theta.original <- theta.original
    v$mplefit <- initialfit
         
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
