############################################################################
# The <ergm.stepping> function provides one of the styles of maximum
# likelihood estimation that can be used. This one is attributed to ?? and
# uses ?? approach. The other  MLE styles are found in functions <ergm.robmon>
# <ergm.stocapprox> and <ergm.mainfitloop>
#
# --PARAMETERS--
#   theta0         : the initial theta values
#   nw             : the network
#   model          : the model, as returned by <ergm.getmodel>
#   Clist          : a list of several network and model parameters,
#                    as returned by <ergm.Cprepare>
#   initialfit     : an ergm object, as the initial fit
#   MCMCparams     : a list of parameters for controlling the MCMC sampling
#   MHproposal     : an MHproposal object for 'nw', as returned by
#                    <getMHproposal>
#   MHproposal.miss: an MHproposal object for the missing network of'nw',
#                    as returned by <getMHproposal>
#   verbose        : whether the MCMC sampling should be verbose AND
#                    the diagnostic plots should be printed ; default=FALSE
#   ...            : additional paramters that are passed onto
#                    <ergm.estimate> and <simulate.formula>
#
# --RETURNED--
#   v: an ergm object as a list containing several items; for details see
#      the return list in the <ergm> function header (<ergm.stepping>=@)
#
###########################################################################      

ergm.stepping = function(theta0, nw, model, Clist, initialfit, 
                         MCMCparams, MHproposal, MHproposal.miss, 
                         verbose=FALSE, ...){

  #   preliminary, to set up structure. 
  nw.orig <- nw
  asyse=theta0-theta0
  mc.se=1+0.05*asyse
  mle.lik=initialfit$mle.lik
  theta.original=theta0
  
  ## Prepare the output structure:
  formula <- model$formula  # formula for this model
  obsstats <- summary(model$formula)  # Observed statistics
  theta0 <- theta0  # beginning parameter value
  samples <- list()  # matrices of sampled network statistics
  sampmeans <- list() # vectors of column means of stats matrices
  xi <- list() # "new obsstats" values, somewhere between obsstats and sampmeans
  eta <- list() # MLE using lognormal approximation and "new obsstats"
  gamma <- list() # factor controlling convex combo: # xi=gamma*obsstats + (1-gamma)*sampmeans	
	
	iter <- 0
	eta[[1]] <- theta0
	finished <- FALSE
  countdown <- 2
	while (!finished) { # Iterate until gamma==1
		iter=iter+1
    ## Generate an mcmc sample from the probability distribution determined by orig.mle
		samples[[iter]]=simulate.formula(formula, nsim=MCMCparams$stepMCMCsize,
                                     theta0=eta[[iter]], burnin=MCMCparams$burnin, 
                                     interval=MCMCparams$interval, statsonly=TRUE,
                                     ...)
		sampmeans[[iter]]=colMeans(samples[[iter]])
		
		hi <- MCMCparams$gridsize  # Goal: Let gamma be largest possible multiple of .01
		lo <- gamm <- 0
		cat("Iteration #",iter, ". ")
    if (verbose) {
      cat("Current canonical parameter:\n")
      print(eta[[iter]])
    }
		while (hi-lo>1 || hi > gamm) {
      gamm<-ceiling((hi+lo)/2)
			gamma[[iter]] = gamm/MCMCparams$gridsize
			xi[[iter]] = gamma[[iter]]*obsstats  + (1-gamma[[iter]])*sampmeans[[iter]]
			inCH=is.inCH(xi[[iter]], samples[[iter]])
			if (inCH) {
        lo <-gamm
        # If 1/gridsize is not small enough, function gives warning message
      } else { 
        hi <- gamm
        if (gamm==1) {
          warning("gamma=", 1/MCMCparams$gridsize, " still not small enough to ",
                  "stay in convex hull.  A larger gridsize than ",
                  MCMCparams$gridsize, " might help.")
        }
      }
		}
		if (!inCH && gamm>1) {# Last attempt was a fail, so decrease gamm by one
			gamm <- gamm-1
			gamma[[iter]] <- gamm/MCMCparams$gridsize
			xi[[iter]] <- gamma[[iter]]*obsstats  + (1-gamma[[iter]])*sampmeans[[iter]]
		}
	
    # Now we have found a gamm that moves xi inside the convex hull, but very 
    # close to the boundary.
    # As in Hummel, Hunter, Handcock, we replace gamm by C*gamm for some 
    # constant C (the paper uses 0.95) and update the xi[[iter]] value prior 
    # to use, for the MCMC MLE step, a xi that we are certain is not on the 
    # edge of the convex hull (except when gamma=1):
    if (gamm == MCMCparams$gridsize && 
        is.inCH((1/.95)*obsstats + (1-1/.95)*sampmeans[[iter]], samples[[iter]]) ) {
      if (verbose) 
        cat("Observed stats are well inside the convex hull.\n")
      C <- 1
      countdown <- countdown - 1
    } else {
      C <- .95
      countdown <- 2 #  Reset countdown
    }
    # We'd like to have gamma==1 for 2 consecutive iterations before
    # we declare that we're finished.
    finished = (countdown==0)
    
    xi[[iter]] <- C*gamma[[iter]]*obsstats  + (1-C*gamma[[iter]])*sampmeans[[iter]]
    # When the stepped xi is in the convex hull (but not on the boundary), find the MLE for gyobs=xi
    #		cat("  Trying gamma=", gamm,"/", MCMCparams$gridsize,"\n") #Only if C=1
		cat("  Trying gamma=", round(C*gamma[[iter]],3),"\n")  
		flush.console()
    
############# PLOTS print if VERBOSE=TRUE #################
  if(verbose){
#  Take a look at obsstats (in red dot) and the "new obsstats" (in green triangle):
##		par(mgp = c(2,.8,0), mar = .1+c(3,3,3,1)) ## How do I put more margin at the top?
    par(ask=TRUE)
		pairs(rbind(samples[[iter]], obsstats, xi[[iter]], sampmeans[[iter]]), 
			  col=c(rep(1, nrow(samples[[iter]])), 2, 2, 2),#2, 7, 3), 
			  pch=c(rep(46, nrow(samples[[iter]])), 4, 20, 4),
#			  main=paste("Iteration ", iter, ":", "\n",
#						 "Stepped stats = ", round(C*gamma[[iter]],3), 
#						 "(Observed stats) + ", round(1-C*gamma[[iter]],3), "(Sample Mean)", sep=""),
			  main=paste("Iteration ", iter, ": ",
                   "gamma = ", round(C*gamma[[iter]], 3), sep=""),
			  cex.main=1.5, cex.axis=1.5, cex.lab=1.5)
    
#		if (!is.null(filename.prefix)) {
#			pdf(paste(filename.prefix, ".iter",iter,".pdf",sep=""))
#			pairs(rbind(samples[[iter]], obsstats, xi[[iter]], sampmeans[[iter]]), 
#				  col=c(rep(1, nrow(samples[[iter]])), 2, 7, 3), 
#				  pch=c(rep(20, nrow(samples[[iter]])),8 ,8 ,8 ),
#				  main=paste("Iteration ", iter, ":  Yellow = ", round(gamma[[iter]],3), 
#							 "(Red) + ", round(1-gamma[[iter]],3), "(Green)", sep=""),
#				  cex.main=1.5, cex.axis=1.5, cex.lab=1.5)
#			dev.off()    
#		}
		} #ends if(verbose)
############# END PLOTS #################		
		
    # ergm.estimate requires that the simulated stats be "centered" in the sense that the
    # observed statistics would be exactly zero on the same scale.  In this case, the
    # "observed statistics" equal xi[[iter]].
    v<-ergm.estimate(theta0=eta[[iter]], model=model, 
                     statsmatrix=sweep(samples[[iter]], 2, xi[[iter]], '-'), 
                     nr.maxit=MCMCparams$nr.maxit,
                     metric=MCMCparams$metric,
                     verbose=verbose,
                     trace=0,  # suppress 'optim' output
                     estimateonly=TRUE, 
                     #statsmatrix.miss=statsmatrix.miss, 
                     #epsilon=MCMCparams$epsilon,
                     # nr.reltol=MCMCparams$nr.reltol,
                     #calc.mcmc.se=MCMCparams$calc.mcmc.se, hessianflag=MCMCparams$hessian,
                     # trustregion=MCMCparams$trustregion, method=MCMCparams$method, 
                     #compress=MCMCparams$compress, 
                     ...)
    eta[[iter+1]]<-v$coef
	}
	cat("Now ending with one large sample for MLE. \n")
	flush.console()
	iter=iter+1
    finalsample=simulate.formula(formula, nsim=MCMCparams$samplesize,
                                     theta0=eta[[iter]], burnin=MCMCparams$burnin, 
                                     interval=MCMCparams$interval, statsonly=TRUE, 
                                     ...)
	sampmeans[[iter]]=colMeans(finalsample)
  # Using the large sample from xi=.95*obsstats+.05*sampmeans we now find the MLE for the original problem:		
	xi[[iter]] = obsstats
	v<-ergm.estimate(theta0=eta[[iter]], model=model, 
									 statsmatrix=sweep(finalsample, 2, xi[[iter]], '-'), 
									 nr.maxit=MCMCparams$nr.maxit,
									 metric=MCMCparams$metric,
									 verbose=verbose,
                   trace=0,  # suppress 'optim' output
									 #estimateonly=TRUE,
                   #statsmatrix.miss=statsmatrix.miss, 
                   epsilon=MCMCparams$epsilon,
                    nr.reltol=MCMCparams$nr.reltol,
                   calc.mcmc.se=MCMCparams$calc.mcmc.se, hessianflag=MCMCparams$hessian,
                    trustregion=MCMCparams$trustregion, method=MCMCparams$method, 
                   compress=MCMCparams$compress, 
									 ...)
	eta[[iter+1]]<-v$coef 
	
############# FINAL PLOT 1 prints if VERBOSE=TRUE #################
	if(verbose){		
	pairs(rbind(finalsample, obsstats, sampmeans[[iter]]), 
		  col=c(rep(1, nrow(finalsample)), 7, 3), 
		  pch=c(rep(46, nrow(finalsample)),4 ,4 ),
		  main=paste("Final Stepping Iteration (#", iter, ")"),# "\n",
#					 "Green = mean, Yellow=observed", sep=""),
          cex.main=1.5, cex.axis=1.5, cex.lab=1.5)
	} #ends if(verbose)
############## END PLOT #################		
#					
#  # Here is the second "larger sample size" loop, using the new (true) MLE as eta0.  The idea is
#  # to both get a slightly better MLE and to get a much better Hessian.
#
#	iter=iter+1
#	finalsample2=simulate.formula(formula, nsim=MCMCparams$samplesize,
#									 theta0=eta[[iter]], burnin=MCMCparams$burnin, 
#									 interval=MCMCparams$interval, statsonly=TRUE, ...)
#					
#	sampmeans[[iter]]=colMeans(finalsample2)
#
############## FINAL PLOT 2 prints if VERBOSE=TRUE #################
#		if(verbose){		
#	pairs(rbind(finalsample2, obsstats, sampmeans[[iter]]), 
#		  col=c(rep(1, nrow(finalsample2)), 7, 3), 
#		  pch=c(rep(46, nrow(finalsample2)),4 ,4 ),
#		  main=paste("Final Estimation"),#:  Green = mean, Yellow=observed", sep=""),
#          cex.main=1.5, cex.axis=1.5, cex.lab=1.5)
#		} #ends if(verbose)
############## END PLOT #################		
#
#  inCH=is.inCH(obsstats, finalsample2)
#	if (!inCH) {
#    stop("Observed statistics are not in final sample convex hull.")
#  }
#
#  # ergm.estimate requires that the simulated stats be "centered" in the sense that the
#  # observed statistics would be exactly zero on the same scale.  In this case, the
#  # "observed statistics" equal obsstats.
#	v<-ergm.estimate(theta0=eta[[iter]], model=model, 
#                   statsmatrix=sweep(finalsample2, 2, obsstats, '-'),
#                   #statsmatrix.miss=statsmatrix.miss, 
#                   epsilon=MCMCparams$epsilon,
#                   nr.maxit=MCMCparams$nr.maxit,
#                   nr.reltol=MCMCparams$nr.reltol,
#                   calc.mcmc.se=MCMCparams$calc.mcmc.se, hessianflag=MCMCparams$hessian,
#                   trustregion=MCMCparams$trustregion, method=MCMCparams$method, 
#                   metric=MCMCparams$metric,
#                   compress=MCMCparams$compress, 
#                   verbose=verbose, 
#                   trace=0,  # suppress 'optim' output
#                   ...)

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


