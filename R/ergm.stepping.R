ergm.stepping = function(theta0, nw, model, Clist, initialfit, 
				MCMCparams=MCMCparams, 
				MHproposal=MHproposal, MHproposal.miss=MHproposal.miss, 
				verbose=FALSE, estimate=TRUE, sequential=TRUE, ...)
{ 
#   preliminary, to set up structure. 
					nw.orig <- nw
					asyse=theta0-theta0
					mc.se=1+0.05*asyse
					mle.lik=initialfit$mle.lik
					theta.original=theta0
##Trash these?	
#					v=list(coef=theta0)
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
										 interval=MCMCparams$interval, statsonly=TRUE)
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
    # If 1/gridsize is not small enough, function fails and gives error message
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
    # When the stepped xi is in the convex hull, find the MLE for gyobs=xi
		cat("  Trying alpha=", alph,"/", MCMCparams$gridsize,"\n")
    flush.console()

    # ergm.estimate requires that the simulated stats be "centered" in the sense that the
    # observed statistics would be exactly zero on the same scale.  In this case, the
    # "observed statistics" equal xi[[iter]].
    
    # As in Hummel, Hunter, Handcock, we replace alph by C*alph for some constant C (the 
    # paper uses 0.95) and update the xi[[iter]] value prior to use, for the MCMC MLE step, 
    # a xi that we are certain is not on the edge of the convex hull:
	C=.95
	xi[[iter]] = C*alpha[[iter]]*obsstats  + (1-C*alpha[[iter]])*sampmeans[[iter]]
		
    v<-ergm.estimate(theta0=eta[[iter]], model=model, 
                     statsmatrix=sweep(samples[[iter]], 2, xi[[iter]], '-'), 
                     nr.maxit=MCMCparams$nr.maxit,
                     metric=MCMCparams$metric,
                     verbose=verbose,
                     estimateonly=TRUE, 
					 #statsmatrix.miss=statsmatrix.miss, 
					 #epsilon=MCMCparams$epsilon,
					 # nr.reltol=MCMCparams$nr.reltol,
					 #calc.mcmc.se=MCMCparams$calc.mcmc.se, hessianflag=MCMCparams$hessian,
					 # trustregion=MCMCparams$trustregion, method=MCMCparams$method, 
					 #compress=MCMCparams$compress, 
                     ...)
    eta[[iter+1]]<-v$coef
		
## If alpha is still not 1, go back to next iteration
## If alpha is 1, then we might be able to obtain a valid MLE using C=1. 
## So we do one more loop with a larger sample size for a xi based on C=.95*alpha=1, 
## (to obtain a better sample),		
## followed by one more loop with a larger sample size for the true xi (C=1, alpha=1)
## (to obtain a better MLE):
	    finished = (alph==MCMCparams$gridsize)
	}
	cat("Now ending with one large sample for MLE. \n")
	flush.console()
	iter=iter+1
    finalsample=simulate.formula(formula, nsim=MCMCparams$samplesize,
                                     theta0=eta[[iter]], burnin=MCMCparams$burnin, 
                                     interval=MCMCparams$interval, statsonly=TRUE)
  # Using the large sample from xi=.95*obsstats+.05*sampmeans we now find the MLE for the original problem:		
	xi[[iter]] = obsstats
	v<-ergm.estimate(theta0=eta[[iter]], model=model, 
									 statsmatrix=sweep(finalsample, 2, xi[[iter]], '-'), 
									 nr.maxit=MCMCparams$nr.maxit,
									 metric=MCMCparams$metric,
									 verbose=verbose,
									 estimateonly=TRUE, 
                                     #statsmatrix.miss=statsmatrix.miss, 
                                     #epsilon=MCMCparams$epsilon,
                                     # nr.reltol=MCMCparams$nr.reltol,
                                     #calc.mcmc.se=MCMCparams$calc.mcmc.se, hessianflag=MCMCparams$hessian,
                                     # trustregion=MCMCparams$trustregion, method=MCMCparams$method, 
                                     #compress=MCMCparams$compress, 
									 ...)
	eta[[iter+1]]<-v$coef 
					
  # Here is the second "larger sample size" loop, using the new (true) MLE as eta0.  The idea is
  # to both get a slightly better MLE and to get a much better Hessian.

	iter=iter+1
	finalsample2=simulate.formula(formula, nsim=MCMCparams$samplesize,
									 theta0=eta[[iter]], burnin=MCMCparams$burnin, 
									 interval=MCMCparams$interval, statsonly=TRUE)			
					
	sampmeans[[iter]]=colMeans(finalsample2)
	inCH=is.inCH(obsstats, finalsample2)
	if (!inCH)
    stop("Observed statistics are not in final sample convex hull.")
    

    # ergm.estimate requires that the simulated stats be "centered" in the sense that the
    # observed statistics would be exactly zero on the same scale.  In this case, the
    # "observed statistics" equal obsstats.
	v<-ergm.estimate(theta0=eta[[iter]], model=model, 
		   statsmatrix=sweep(finalsample2, 2, obsstats, '-'),
           #statsmatrix.miss=statsmatrix.miss, 
           epsilon=MCMCparams$epsilon,
           nr.maxit=MCMCparams$nr.maxit,
           nr.reltol=MCMCparams$nr.reltol,
           calc.mcmc.se=MCMCparams$calc.mcmc.se, hessianflag=MCMCparams$hessian,
           trustregion=MCMCparams$trustregion, method=MCMCparams$method, 
           metric=MCMCparams$metric,
           compress=MCMCparams$compress, 
           verbose=verbose, ...)

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
	v$null.deviance <- 2*network.dyadcount(nw.orig)*log(2)
	v$mle.lik <- mle.lik
	v$etamap <- model$etamap

	v
}  # Ends the whole function


