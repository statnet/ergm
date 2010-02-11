ergm.stepping = function(theta0, nw, model, Clist, initialfit, 
#control=control.ergm(nsim1=100, nsim2=1000, gridsize=100),  # simulation parameters
## approx=lognormapprox, filename.prefix=NULL, plots=FALSE,  # currently useless, but plots can be reimplemented
				MCMCparams=MCMCparams, 
				MHproposal=MHproposal, MHproposal.miss=MHproposal.miss, 
				verbose=FALSE, estimate=TRUE, sequential=TRUE
#,...
){ 

#   preliminary, to set up structure. 
					nw.orig <- nw
					asyse=theta0-theta0
					mc.se=1+0.05*asyse
					mle.lik=initialfit$mle.lik
#					v=list(coef=theta0)
					theta.original=theta0
#					thetaprior=theta0-theta0
#					stats <- matrix(0,ncol=Clist$nstats,nrow=MCMCparams$samplesize)					
#					stats[1,] <- Clist$obs - Clist$meanstats
#					MCMCparams$stats <- stats
#					MCMCparams$meanstats <- Clist$meanstats
#					thetaprior <- theta0
#					theta0 <- v$coef
#					eta0 <- ergm.eta(theta0, model$etamap)
						
	
	
## Prepare the output structure:
   	formula = model$formula  # formula for this model
    obsstats= summary(model$formula)  # Observed statistics
    theta0 = theta0  # beginning parameter value
    samples=list()  # matrices of sampled network statistics
    sampmeans=list() # vectors of column means of stats matrices
    xi=list() # "new obsstats" values, somewhere between obsstats and sampmeans
    eta=list() # MLE using lognormal approximation and "new obsstats"
    alpha=list() # factor controlling convex combo: # xi=alpha*obsstats + (1-alpha)*sampmeans	
	
	iter = 0
	eta[[1]]=theta0
	finished = FALSE
	while (!finished) { # Iterate until alpha==1
		iter=iter+1
## Generate an mcmc sample from the probability distribution determined by orig.mle
		samples[[iter]]=simulate.formula(formula, nsim=MCMCparams$stepMCMCsize,       
										 theta0=eta[[iter]], burnin=MCMCparams$burnin, 
										 interval=MCMCparams$interval, statsonly=T)
		sampmeans[[iter]]=colMeans(samples[[iter]])
		
		hi <- MCMCparams$gridsize  # Goal: Let alpha be largest possible multiple of .01
		lo <- alph <- 0
		cat("Iteration #",iter, ". ")
		while (hi-lo>1 || hi > alph) {
      alph<-ceiling((hi+lo)/2)
			alpha[[iter]] = alph/MCMCparams$gridsize
			xi[[iter]] = alpha[[iter]]*obsstats  + (1-alpha[[iter]])*sampmeans[[iter]]
			inCH=is.inCH(xi[[iter]], samples[[iter]])
			if (inCH)
			lo=alph
## If 1/gridsize is not small enough, function fails and gives error message
			else if (alph==1)
			stop("alpha=", 1/MCMCparams$gridsize, "still not small enough for MLE.N.  Try using a", 
				 " gridsize larger than ", MCMCparams$gridsize)
			else
			hi=alph
		}
		if (!inCH) {# Last attempt was a fail, so decrease alph by one
			alph = alph-1
			alpha[[iter]] = alph/MCMCparams$gridsize
			xi[[iter]] = alpha[[iter]]*obsstats  + (1-alpha[[iter]])*sampmeans[[iter]]
		}
## When the stepped xi is in the convex hull, find the MLE for gyobs=xi
##########Here is the problem:  why do I not get the same estimate by each method? ##########
####		cat(paste("xi",xi[[iter]], "\n"))
####    eta[[iter+1]]=optim(par=eta[[iter]], approx, control=list(fnscale=-1), eta0=eta[[iter]], s=samples[[iter]], xi=xi[[iter]])$par
####		cat(paste("eta", eta[[iter+1]], "\n"))
		cat("  Trying alpha=", alph,"/", MCMCparams$gridsize,"\n")
    flush.console()

    v<-ergm.estimate(theta0=eta[[iter]], model=model, 
                     xobs=xi[[iter]] - sampmeans[[iter]],  
						 statsmatrix=samples[[iter]], 
#statsmatrix.miss=statsmatrix.miss, 
#epsilon=MCMCparams$epsilon,
#nr.maxit=MCMCparams$nr.maxit,
# nr.reltol=MCMCparams$nr.reltol,
#calc.mcmc.se=MCMCparams$calc.mcmc.se, hessian=MCMCparams$hessian,
# trustregion=MCMCparams$trustregion, method=MCMCparams$method, metric="Likelihood",
#compress=MCMCparams$compress, 
						 verbose=verbose,
						 estimateonly=TRUE)		
    eta[[iter+1]]<-v$coef
		
## If alpha is still not 1, go back to next iteration
## If alpha is 1, then latest eta is a valid MLE, but maybe one more loop
## would not hurt with a larger sample size
	    finished = (alph==MCMCparams$gridsize)
	}
	cat("Now ending with one large sample for MLE. \n")
	flush.console()
	iter=iter+1
    finalsample=simulate.formula(formula, nsim=MCMCparams$samplesize,   ###why can I not change this to MCMCparams$MCMCsamplesize without an error in length?
                                     theta0=eta[[iter]], burnin=MCMCparams$burnin, 
                                     interval=MCMCparams$interval, statsonly=T)
	sampmeans[[iter]]=colMeans(finalsample)
	inCH=is.inCH(obsstats, finalsample)
	if (!inCH)
    stop("Observed statistics are not in final sample convex hull.")
    
## final iteration now that we are in the convex hull, needs to be same form as what is returned in ergm.estimate
#####	final.mle = eta[[iter+1]] = optim(par=eta[[iter]], approx,
#####									  control=list(fnscale=-1), eta0=eta[[iter]], 
#####									  s=finalsample, xi=obsstats)$par

	v<-ergm.estimate(theta0=eta[[iter]], model=model, 
                   xobs=obsstats - sampmeans[[iter]],  
					 statsmatrix=finalsample,
#statsmatrix.miss=statsmatrix.miss, 
#epsilon=MCMCparams$epsilon,
#nr.maxit=MCMCparams$nr.maxit,
#nr.reltol=MCMCparams$nr.reltol,
#calc.mcmc.se=MCMCparams$calc.mcmc.se, hessian=MCMCparams$hessian,
#trustregion=MCMCparams$trustregion, method=MCMCparams$method, metric="Likelihood",
#compress=MCMCparams$compress, 
									 verbose=verbose)
	
#####	final.mle
	mle.lik <- mle.lik + abs(v$loglikelihood)
	v$newnetwork <- nw
	v$burnin <- MCMCparams$burnin
	v$samplesize <- MCMCparams$samplesize
	v$interval <- MCMCparams$interval
	v$network <- nw.orig
	v$newnetwork <- nw
	v$interval <- MCMCparams$interval
	v$theta.original <- theta.original
	v$mplefit <- initialfit
	v$parallel <- MCMCparams$parallel
	
	if(!v$failure & !any(is.na(v$coef))){
#     asyse <- sqrt(diag(robust.inverse(-v$hessian)))
#     asyse <- try(sqrt(diag(robust.inverse(cov(statsmatrix)))))
		asyse <- mc.se
		options(warn=-1)
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
# }  #I am really not sure about this mle.lik value.
	v
}  # Ends the whole function


