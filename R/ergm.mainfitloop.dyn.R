ergm.mainfitloop.dyn <- function(theta0, nw, model.form, model.diss, 
                                 Clist,
                                 gamma, initialfit, 
                                 MCMCparams, MHproposal.form, MHproposal.diss,
                                 verbose=FALSE,
                                 epsilon=1e-10,
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

  if(is.null(Clist$meanstats)){Clist$meanstats <- Clist$obs}
  stats <- matrix(0,ncol=Clist$nparam,nrow=MCMCparams$samplesize)
  stats[1,] <- summary(model.form$formula, basis=nw) - Clist$meanstats

  MCMCparams=c(MCMCparams, list(stats=stats,
               meanstats=Clist$meanstats, orig.obs=Clist$obs))

#  while(any(mcmc.precision*asyse < mc.se, na.rm=TRUE) && iteration <= MCMCparams$maxit){
  while(iteration <= MCMCparams$maxit){
    thetaprior <- theta0
    theta0 <- v$coef
    eta0 <- ergm.eta(theta0, model.form$etamap)
    cat("Iteration ", iteration,": Sampling ", MCMCparams$samplesize,
        " with parameter: \n", sep="")
    print(theta0)
#   MCMCparams=list(samplesize=MCMCparams$samplesize, burnin=burnin, interval=interval,parallel=parallel,
#                   gamma=gamma, dyninterval=dyninterval,
#                   meanstats=Clist$meanstats, orig.obs=Clist$obs,
#                   maxchanges=10*maxchanges)
#   MHproposal.form=list(package=proposalpackage, type=proposaltype)
    z <- ergm.getMCMCDynsample(nw, model.form, model.diss,
                               MHproposal.form, MHproposal.diss,
                               eta0, gamma, MCMCparams, verbose)
    statsmatrix <- z$statsmatrix
#   print( summary(nw) )
    nw <- z$newnetwork
# MCMCparams$orig.obs <- summary(model.form$formula, basis=newnetwork)
#   MCMCparams$orig.obs <- statsmatrix[nrow(statsmatrix),]+meanstats
#   Next line VIP for sequential update
    MCMCparams$stats[1,] <- summary(model.form$formula, basis=nw)-MCMCparams$meanstats
#   MCMCparams$stats[1,] <- statsmatrix[nrow(statsmatrix),]
#   print( meanstats )
#   print( MCMCparams$orig.obs )
#   print( apply(statsmatrix,2,mean) )
#   print( statsmatrix[nrow(statsmatrix),] )
#   print(  MCMCparams$stats[1,] )
#
#   Clist <- ergm.Cprepare(newnetwork, model.form)
#   Clist$meanstats=meanstats
#   aaa <<- nw
#   print(summary(newnetwork ~ hamming(aaa)))
#   aaa <<- g0
#   print(summary(newnetwork ~ hamming(aaa)))
#
    if(verbose){cat("Back from unconstrained MCMC...\n")}
    iteration <- iteration + 1
    if(MCMCparams$steplength<1 && iteration < MCMCparams$maxit ){
      statsmean <- apply(statsmatrix,2,mean)
      statsmatrix <- sweep(statsmatrix,2,(1-MCMCparams$steplength)*statsmean,"-")
    }
    if(ergm.checkdegeneracy(statsmatrix, verbose=verbose)){
     if(iteration <= MCMCparams$maxit){
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
                network.density(nw),"\n"))
      cat(paste("The density of the original network is",
                network.density(nw.orig),"\n"))
      cat("Summary of the simulation:\n")
      print(apply(statsmatrix,2,summary.statsmatrix.ergm),scipen=6)
      degreedist(nw)
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
              samplesize=MCMCparams$samplesize, failure=TRUE,
              newnetwork = nw,
              formula = model.form$formula)
    return(structure (l, class="ergm"))
  }
  if(verbose){cat("Calling MCMLE Optimization...\n")}
  if(verbose){cat("Using Newton-Raphson Step ...\n")}
#
# If not the last iteration do not compute all the extraneous
# statistics that are not needed until output
#
  if(iteration <= MCMCparams$maxit){
   v<-ergm.estimate.only(theta0=theta0, model=model.form,
                    statsmatrix=statsmatrix,
                    statsmatrix.miss=NULL,
                    epsilon=epsilon,
                    nr.maxit=MCMCparams$nr.maxit, 
                    calc.mcmc.se=MCMCparams$calc.mcmc.se, 
                    hessian=MCMCparams$hessian,
                    trustregion=MCMCparams$trustregion, 
                    method=MCMCparams$method, metric="Likelihood",
                    compress=MCMCparams$compress, verbose=verbose)
  }
#
# End main loop
#
  }
# This is the last iteration, so compute all the extraneous
# statistics that are not needed until output
#
  v<-ergm.estimate(theta0=theta0, model=model.form,
                   statsmatrix=statsmatrix, 
                   statsmatrix.miss=NULL,
                   epsilon=epsilon,
                   nr.maxit=MCMCparams$nr.maxit, 
                   calc.mcmc.se=MCMCparams$calc.mcmc.se, 
                   hessian=MCMCparams$hessian,
                   trustregion=MCMCparams$trustregion, 
                   method=MCMCparams$method, metric="Likelihood",
                   compress=MCMCparams$compress, verbose=verbose)
#
#    Reform the output
#
# check for degeneracy?  See ergm.checkdegeneracy
# l$boundDeg <- list(condAllDegExact=condAllDegExact, maxout=maxout,
#                    maxin=maxin, minout=minout,
#                    minin=minin, attribs=attribs)
#   This is the MCMC estimate as in the Hunter&Handcock paper
#  l$mcmcloglik <- l$mcmcloglik - Clist$ndyads*log(2)
#   This is the MPLE value estimate plus the ratio
#   set in the last line of "ergm.R" 
#   v$mle.lik <- -glm.fit$deviance/2 + v$loglikelihood
  mle.lik <- mle.lik + v$loglikelihood

  v$newnetwork <- nw
  v$formula <- model.form$formula

###### old ergm.statseval ends here    
    
    v$burnin <- MCMCparams$burnin
    v$samplesize <- MCMCparams$samplesize
    v$interval <- MCMCparams$interval
    v$network <- nw.orig
    v$newnetwork <- nw
    v$interval <- MCMCparams$interval
    v$theta.original <- theta.original
    v$proposal.form <- MHproposal.form$name
    v$proposal.diss <- MHproposal.diss$name
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
  v$mcmcloglik <- v$mcmcloglik - network.dyadcount(nw.orig)*log(2)
# if(!is.na(v$mle.lik) && v$mle.lik>0){
#   v$mle.lik <- -v$mplefit$glm$deviance/2 
# }  #I'm really not sure about this mle.lik value.
  v
}
