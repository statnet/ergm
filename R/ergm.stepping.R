#  File R/ergm.stepping.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2015 Statnet Commons
#######################################################################
############################################################################
# The <ergm.stepping> function provides one of the styles of maximum
# likelihood estimation that can be used. This one is attributed to ?? and
# uses ?? approach. The other  MLE styles are found in functions <ergm.robmon>
# <ergm.stocapprox> and <ergm.mainfitloop>
#
# --PARAMETERS--
#   init         : the initial theta values
#   nw             : the network
#   model          : the model, as returned by <ergm.getmodel>
#   Clist          : a list of several network and model parameters,
#                    as returned by <ergm.Cprepare>
#   initialfit     : an ergm object, as the initial fit
#   control     : a list of parameters for controlling the MCMC sampling
#   MHproposal     : an MHproposal object for 'nw', as returned by
#                    <MHproposal>
#   MHproposal.obs: an MHproposal object for the observed network of'nw',
#                    as returned by <MHproposal>
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

ergm.stepping = function(init, nw, model, initialfit, constraints,
                         control, MHproposal, MHproposal.obs, 
                         verbose=FALSE, ...){

  #   preliminary, to set up structure. 
  nw.orig <- nw
  asyse=init-init
  mc.se=1+0.05*asyse
  mle.lik=initialfit$mle.lik
  theta.original=init
  
  ## Prepare the output structure:
  formula <- model$formula  # formula for this model
  obsstats <- summary(model$formula)  # Observed statistics
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
		samples[[iter]]=simulate.formula(formula, nsim=control$Step.MCMC.samplesize,
                                     coef=eta[[iter]], statsonly=TRUE,
                                     constraints=constraints, 
                                     control=set.control.class("control.simulate.formula",control), ...)
		sampmeans[[iter]]=colMeans(samples[[iter]])
		
		hi <- control$Step.gridsize  # Goal: Let gamma be largest possible multiple of .01
		lo <- gamm <- 0
		cat("Iteration #",iter, ". ")
    if (verbose) {
      cat("Current canonical parameter:\n")
      print(eta[[iter]])
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
        cat("Observed stats are well inside the convex hull.\n")
      countdown <- countdown - 1
    } else {
      countdown <- 2 #  Reset countdown
    }
    # We'd like to have gamma==1 for 2 consecutive iterations before
    # we declare that we're finished.
    finished = (countdown==0) || (iter >= control$Step.maxit)
    
    # When the stepped xi is in the convex hull (but not on the boundary), find the MLE for gyobs=xi
		cat("  Trying gamma=", gamma[[iter]],"\n")  
		flush.console()
    
    ############# PLOTS print if VERBOSE=TRUE #################
    if(verbose){
      # Take a look at obsstats (in red dot) and the "new obsstats" (in green triangle):
      # par(mgp = c(2,.8,0), mar = .1+c(3,3,3,1)) ## How do I put more margin at the top?
      par(ask=TRUE)
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
                     statsmatrix=sweep(samples[[iter]], 2, xi[[iter]], '-'), 
                     nr.maxit=control$MCMLE.NR.maxit,
                     metric=control$MCMLE.metric,
                     verbose=verbose,
                     trace=0,  # suppress 'optim' output
                     estimateonly=TRUE, 
                     #statsmatrix.obs=statsmatrix.obs, 
                     #epsilon=control$epsilon,
                     # nr.reltol=control$MCMLE.NR.reltol,
                     #calc.mcmc.se=control$MCMC.addto.se, hessianflag=control$main.hessian,
                     # trustregion=control$MCMLE.trustregion, method=control$MCMLE.method, 
                     #compress=control$MCMC.compress, 
                     ...)
    eta[[iter+1]]<-v$coef
	}
	cat("Now ending with one large sample for MLE. \n")
	flush.console()
	iter <- iter+1
  finalsample <- simulate.formula(formula, nsim=control$MCMC.samplesize,
                                  coef=eta[[iter]], statsonly=TRUE, 
                                  constraints=constraints, 
                                  control=set.control.class("control.simulate.formula",control), ...)
  sampmeans[[iter]] <- colMeans(finalsample)
  xi[[iter]] <- obsstats
	v<-ergm.estimate(init=eta[[iter]], model=model, 
									 statsmatrix=sweep(finalsample, 2, xi[[iter]], '-'), 
									 nr.maxit=control$MCMLE.NR.maxit,
									 metric=control$MCMLE.metric,
									 verbose=verbose,
                   trace=0,  # suppress 'optim' output
									 #estimateonly=TRUE,
                   #statsmatrix.obs=statsmatrix.obs, 
                   epsilon=control$epsilon,
                    nr.reltol=control$MCMLE.NR.reltol,
                   calc.mcmc.se=control$MCMC.addto.se, hessianflag=control$main.hessian,
                    trustregion=control$MCMLE.trustregion, method=control$MCMLE.method, 
                   compress=control$MCMC.compress, 
									 ...)
  eta[[iter+1]] <- v$coef
	
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
  #####	final.mle
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
	
	if(!v$failure & !any(is.na(v$coef))){
		asyse <- mc.se
		if(is.null(v$covar)){
			asyse[names(v$coef)] <- suppressWarnings(sqrt(diag(ginv(-v$hessian))))
		}else{
			asyse[names(v$coef)] <- suppressWarnings(sqrt(diag(v$covar)))
		}
	}
	
  v$sample <- ergm.sample.tomcmc(v$sample, control)
  v$etamap <- model$etamap
  v$iterations <- iter

	v
}  # Ends the whole function

## Given two matrices x1 and x2 with d columns (and any positive
## numbers of rows), find the greatest gamma<=steplength.max s.t., the
## points of x2 shrunk towards the centroid of x1 a factor of
## margin*gamma, are all in the convex hull of x1.

## This is a variant of Hummel et al. (2013)'s steplength algorithm
## also usable for missing data MLE.
.Hummel.steplength <- function(x1, x2=NULL, margin=0.05, steplength.max=1, steplen=1, verbose=FALSE){
  margin <- 1 + margin
  x1 <- rbind(x1); m1 <- colMeans(x1); x1 <- unique(x1)
  if(is.null(x2)) x2 <- rep(0,ncol(x1))
  x2 <- rbind(x2); x2 <- unique(x2)

  if(nrow(x2) > 100){
#  Reduce matrix to those extremes with the num th smallest probabilities
   if(verbose){cat("Using fast and approximate Hummel et al search.")}
   af <- function(x,num=100){
	  num <- min(nrow(x),num)
	  d2 <- sweep(x,2,colMeans(x))
	  sv <- svd(d2)
	  d = d2 %*% sv$v
	  l2 <- rep(0,length(sv$d))
	  ev <- sv$d*sv$d / mean(sv$d*sv$d)
	  l2[ev>0.001] <- 1/ev[ev>0.001]
	  # a are the Mahalanobis distances
	  a = apply(sweep(d*d,2,l2,"*"),1,sum)
          x[order(-a)[1:num],]
   }

   # is.inCH() can handle 1-row M, so this is perfectly fine.
   zerofn <- function(gamma){
	  x=af(t(margin*gamma * t(x2)  + (1-margin*gamma)*m1))
	  is.inCH(x, x1, numcheck=50, verbose=verbose)
   }
   maxit <- 10
  }else{
   if(verbose){cat("Using detailed Hummel et al search.")}
   maxit <- 25
   zerofn <- function(gamma){is.inCH(t(margin*gamma * t(x2)  + (1-margin*gamma)*m1), x1, verbose=verbose)}
  }

  # Modify search to reflect the binary nature of the outcome and the
  # prior belief about the gamma (centered on prior valeu with a s.d. of sd.p)
  sd.p <- 0.1
  low <- 0
  high <- steplength.max
  g <- steplen
  i <- 0
  while(i < maxit & abs(high-low)>0.01){
   z=zerofn(g)
   if(z){
    low <- g
    g <- qnorm( mean(pnorm(c(high,g), mean = steplen, sd = sd.p)), mean = steplen, sd = sd.p)
   }else{
    high <- g
    g <- qnorm( mean(pnorm(c(low, g), mean = steplen, sd = sd.p)), mean = steplen, sd = sd.p)
   }
   i <- i+1
#  out <- c(i,g,low,high,z)
#  names(out) <- c("iters","est","low","high","z")
#  print(out)
   if(verbose) cat(sprintf("iter= %d, est=%f, low=%f, high=%f, test=%d.\n",i,g,low,high,z))
  }
  g
}
