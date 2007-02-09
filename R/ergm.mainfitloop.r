ergm.mainfitloop <- function(theta0, nw, model, Clist,
                             BD, initialfit, burnin,
                             MCMCsamplesize, interval, maxit, proposaltype,
                             proposalpackage, compress=FALSE, verbose=FALSE,
                             mcmc.precision=0.05, 
                             epsilon=1e-10, nr.maxit=100, calc.mcmc.se=TRUE,
                             hessian=FALSE, trustregion=20,
                             steplength=1,
                             metric="Likelihood",
                             method="BFGS", estimate=TRUE, ...) {
  iteration <- 1
  newnetwork <- nw
  smean <- rep(0,Clist$nparam)
#
# meanstats <- Clist$meanstats
#
  asyse=theta0-theta0
  mc.se=1+0.05*asyse
  mle.lik=initialfit$mle.lik
  v=list(coef=theta0)
  theta.original=theta0
  thetaprior=theta0-theta0

#  while(any(mcmc.precision*asyse < mc.se, na.rm=TRUE) && iteration <= maxit){
  while(iteration <= maxit){
    thetaprior <- theta0
    theta0 <- v$coef
    eta0 <- ergm.eta(theta0, model$etamap)
    cat("Iteration ", iteration,": Sampling ", MCMCsamplesize,
        " with parameter: \n", sep="")
    print(theta0)
    MCMCparams=list(samplesize=MCMCsamplesize, burnin=burnin, interval=interval)
    MHproposal=list(package=proposalpackage, type=proposaltype)
    z <- ergm.getMCMCsample(Clist, model, MHproposal, eta0, MCMCparams, verbose, BD)
    statsmatrix=z$statsmatrix
    newnetwork <- network.update(nw, z$newedgelist)
#
#   Clist <- ergm.Cprepare(newnetwork, model)
#   Clist$meanstats=meanstats
#   aaa <<- nw
#   print(summary(newnetwork ~ hamming(aaa)))
#   aaa <<- g0
#   print(summary(newnetwork ~ hamming(aaa)))
#
    if(verbose){cat("Back from unconstrained MCMC...\n")}
    iteration <- iteration + 1
    if(steplength<1 && iteration < maxit ){
      statsmean <- apply(statsmatrix,2,mean)
      statsmatrix <- sweep(statsmatrix,2,(1-steplength)*statsmean,"-")
    }
    if(ergm.checkdegeneracy(statsmatrix, verbose=verbose)){
     if(iteration <= maxit){
      cat(paste("The MCMC sampler is producing degenerate samples.\n",
                "Try starting the algorithm at an alternative model\n",
                "(That is, changing the 'theta0' argument).\n",
                "I am trying something simple...\n"))
#     shrink <- ergm(nw ~ edges)$coef
#     theta0 <- 0.4*theta0
#     theta0["edges"] <- theta0["edges"] + 0.6*shrink
#     v <- list(coef=theta0)
#     v$coef <- theta0
      v$coef <- 0.9*theta0
      next
     }else{
      warning(paste("\n","The MCMC sampler is producing degenerate samples.\n",
                 "Try starting the algorithm at an alternative model\n",
                 "(That is, changing the 'theta0' argument).\n"))
#     v <- list(coef=theta0)
      v$coef <- theta0
      return(v)
     }
    }
    if(verbose){
      cat(paste("The density of the returned network is",
                network.density(newnetwork),"\n"))
      cat(paste("The density of the original network is",
                network.density(nw),"\n"))
      cat("Summary of the simulation:\n")
      print(apply(statsmatrix,2,summary.statsmatrix.ergm),scipen=6)
      degreedist(newnetwork)
    }
    if(verbose){cat("Calling optimization routines...\n")}

###### old statseval begins here.
    
  if(verbose) {cat("statseval:  epsilon value is ",epsilon,"\n")}
  if(!estimate){
    if(verbose){cat("Skipping optimization routines...\n")}
    l <- list(coef=theta0, mc.se=rep(NA,length=length(theta0)),
              sample=statsmatrix, iterations=1, MCMCtheta=theta0,
              loglikelihood=NA, mcmcloglik=NULL, mle.lik=NULL,
              gradient=rep(NA,length=length(theta0)), acf=NULL,
              samplesize=MCMCsamplesize, failure=TRUE,
              newnetwork = newnetwork,
              formula = model$formula)
    return(structure (l, class="ergm"))
  }
  if(verbose){cat("Calling MCMLE Optimization...\n")}
  if(verbose){cat("Using Newton-Raphson Step ...\n")}
#
# If not the last iteration do not compute all the extraneous
# statistics that are not needed until output
#
  if(iteration <= maxit){
   v<-ergm.estimate.only(theta0=theta0, model=model,
                    statsmatrix=statsmatrix, epsilon=epsilon,
                    nr.maxit=nr.maxit, calc.mcmc.se=calc.mcmc.se, hessian=hessian,
                    trustregion=trustregion, method=method, metric="Likelihood",
                    compress=compress, verbose=verbose)
  }
#
# End main loop
#
  }
# This is the last iteration, so compute all the extraneous
# statistics that are not needed until output
#
  v<-ergm.estimate(theta0=theta0, model=model,
                   statsmatrix=statsmatrix, epsilon=epsilon,
                   nr.maxit=nr.maxit, calc.mcmc.se=calc.mcmc.se, hessian=hessian,
                   trustregion=trustregion, method=method, metric="Likelihood",
                   compress=compress, verbose=verbose)
#
#    Reform the output
#
# check for degeneracy?  See ergm.checkdegeneracy
# l$boundDeg <- list(condAllDegExact=condAllDegExact, maxout=maxout,
#                    maxin=maxin, minout=minout,
#                    minin=minin, attribs=attribs)
  v$boundDeg <- BD
#   This is the MCMC estimate as in the Hunter&Handcock paper
#  l$mcmcloglik <- l$mcmcloglik - Clist$ndyads*log(2)
#   This is the MPLE value estimate plus the ratio
#   set in the last line of "ergm.r" 
#   v$mle.lik <- -glm.fit$deviance/2 + v$loglikelihood
  mle.lik <- mle.lik + v$loglikelihood

  v$newnetwork <- newnetwork
  v$formula <- model$formula

###### old ergm.statseval ends here    
    
    v$burnin <- burnin
    v$samplesize <- MCMCsamplesize
    v$interval <- interval
    v$network <- nw
    v$newnetwork <- newnetwork
    v$interval <- interval
    v$theta.original <- theta.original
    v$proposaltype <- proposaltype
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

  endrun <- burnin+interval*(MCMCsamplesize-1)
  attr(v$sample, "mcpar") <- c(burnin+1, endrun, interval)
  attr(v$sample, "class") <- "mcmc"
# v$null.deviance <- v$mplefit$glm$null.deviance
  v$null.deviance <- 2*network.dyadcount(nw)*log(2)
# v$mle.lik <- -v$mplefit$glm$deviance/2 + v$loglikelihood
# v$mle.lik <- -v$null.deviance/2 + v$loglikelihood
# v$mle.lik <- -initialfit$mle.lik/2 + v$loglikelihood
  v$mle.lik <- mle.lik
  v$mcmcloglik <- v$mcmcloglik - network.dyadcount(nw)*log(2)
# if(!is.na(v$mle.lik) && v$mle.lik>0){
#   v$mle.lik <- -v$mplefit$glm$deviance/2 
# }  #I'm really not sure about this mle.lik value.
  v
}
