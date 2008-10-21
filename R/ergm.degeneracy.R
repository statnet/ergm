#  File ergm/R/ergm.degeneracy.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of Washington
#                David R. Hunter, Penn State University
#                Carter T. Butts, University of California - Irvine
#                Steven M. Goodreau, University of Washington
#                Martina Morris, University of Washington
# Copyright 2007 The statnet Development Team
######################################################################
ergm.degeneracy <- function(object, 
                          control=control.ergm(),
                          fast=TRUE,
                          test.only=FALSE,
                          verbose=FALSE) {
  
  if(!is.ergm(object)){
    stop("A ergm object argument must be given.")
  }
  if(is.matrix(object$sample)){
   if(is.null(object$mplefit$glm)){
    current.warn <- options()$warn
    options(warn=-1)
    fit <- try(ergm(object$formula, MPLEonly=TRUE, Mlestimate=FALSE),silent=TRUE)
    options(warn=current.warn)
    if(inherits(fit,"try-error")){
     object$degeneracy.value <- NA
     object$degeneracy.type <- NULL
     return(invisible(object))
    }
    if(is.null(fit$mplefit)){
     object$mplefit$glm <- fit$glm
    }else{
     object$mplefit <- fit$mplefit
    }
   }
   # So a MCMC fit
   if(object$loglikelihood>control$trustregion-0.1){
    object$degeneracy.value <- Inf
   }else{
    changeobs <- (-2*object$mplefit$glm$y+1)*model.matrix(object$mplefit$glm)
    if(fast && nrow(changeobs) > 100){
     index <- sample((1:nrow(changeobs)), size=100, replace=FALSE)
     changeobs <- changeobs[index,]
     wgts <- object$mplefit$glm$prior.weights[index]
    }else{
     wgts <- object$mplefit$glm$prior.weights
     if(nrow(changeobs) > 1000){
      cat("This computation may take a while ...\n")
     }
    }
    object$degeneracy.type <- try(
      apply(changeobs, 1, ergm.compute.degeneracy,
      theta0=object$MCMCtheta, etamap=object$etamap, 
      statsmatrix=object$sample[,!object$etamap$offsettheta,drop=FALSE],
      trustregion=control$trustregion),silent=TRUE)
    if(inherits(object$degeneracy.type,"try-error")){
     object$degeneracy.value <- Inf
     object$degeneracy.type <- NULL
    }else{
     object$degeneracy.type <- t(rbind(object$degeneracy.type,wgts))
     colnames(object$degeneracy.type)[ncol(object$degeneracy.type)] <- "num.dyads"
     object$degeneracy.value <- max(object$degeneracy.type[,1],na.rm=TRUE)
    }
   }
  }else{
   # So a non-MCMC fit
   if("glm" %in% class(object$glm)){
   # So the MPLE was fit
    # This is the change in log-likelihood for logistic regression
#   object$degeneracy.type <- abs(model.matrix(object$glm) %*% object$glm$coef)
    changebeta <- t(influence(object$glm,do.coef=TRUE)$coefficients/object$glm$prior.weights)
#   newbeta <- sweep(changebeta,1,object$glm$coef,"+")
#   changexbeta <- diag(model.matrix(object$glm) %*% newbeta)
    changexchangebeta <- as.matrix(model.matrix(object$glm)) %*% changebeta
    pi <- predict(object$glm,type="response")
    changexpi <- pi %*% as.matrix(model.matrix(object$glm)) 
    changenorm <- as.vector(pi %*% changexchangebeta)*object$glm$prior.weights
    changeobs <- as.vector(object$glm$y %*% changexchangebeta)
    changesum <- changeobs * object$glm$prior.weights
    changey <- as.vector((2*object$glm$y-1) * diag(changexchangebeta))
#   changeobs <- changexbeta %*% (object$glm$prior.weights*object$glm$y)
#   object$degeneracy.type <- abs(sum(changeobs) - changexbeta*(2*object$glm$y-1))
    object$degeneracy.type <- abs(changesum-changey-changenorm)
#   object$degeneracy.type <- changey
    wgts <- object$glm$prior.weights
    object$degeneracy.type <- cbind(object$degeneracy.type,wgts)
    colnames(object$degeneracy.type) <- c("delta.log.lik","num.dyads")
    object$degeneracy.value <- max(object$degeneracy.type[,1],na.rm=TRUE)
   }else{
    object$degeneracy.value <- Inf
    object$degeneracy.type <- NULL
   }
  }
  if(object$degeneracy.value>control$trustregion-0.1){
   object$degeneracy.value <- Inf
  }
  if(is.infinite(object$degeneracy.value)){
   cat("\n Warning: The diagnostics indicate that the model is very unstable.\n   They suggest that the model is degenerate,\n   and that the numerical summaries are suspect.\n")
  }else{
    if(!test.only || object$degeneracy.value > 1){
     cat("The instability of the model is: ",
        format(object$degeneracy.value, digits=2),"\n")
    }
    if(object$degeneracy.value > 1){
      cat("Instabilities greater than 1 suggest the model is degenerate.\n")
    }
  }
  if(verbose){
    print(object$degeneracy.type)
  }
  return(invisible(object))
}

ergm.compute.degeneracy<-function(xobs, theta0, etamap, statsmatrix,
                        epsilon=1e-10, nr.maxit=100, nr.reltol=0.01,
                        verbose=FALSE, trace=6*verbose,
                        hessian=FALSE, guess=theta0,
                        trustregion=20, ...) {
  samplesize <- dim(statsmatrix)[1]
  statsmatrix0 <- statsmatrix
  probs <- rep(1/nrow(statsmatrix0),nrow(statsmatrix0))
  statsmatrix0.miss <- NULL
  probs.miss <- NULL
  av <- apply(sweep(statsmatrix0,1,probs,"*"), 2, sum)
  xsim <- sweep(statsmatrix0, 2, av,"-")
  xsim.miss <- NULL
  probs.miss <- NULL
# xobs0 <- summary(model$formula)
  xobs <- -xobs - av
#
# Set up the initial estimate
#
  if (verbose) cat("Converting theta0 to eta0\n")
  eta0 <- ergm.eta(theta0, etamap) #unsure about this
  etamap$theta0 <- theta0
#
# Log-Likelihood and gradient functions
#
  penalty <- 0.5
  if (verbose) cat("Optimizing loglikelihood\n")
  Lout <- try(optim(par=guess, 
                    fn=llik.fun, #  gr=llik.grad,
                    hessian=FALSE,
                    method="BFGS",
                    control=list(trace=trace,fnscale=-1,reltol=nr.reltol,
                                 maxit=nr.maxit),
                    xobs=xobs,
                    xsim=xsim, probs=probs,
                    xsim.miss=xsim.miss, probs.miss=probs.miss,
                    penalty=0.5, trustregion=trustregion,
                    eta0=eta0, etamap=etamap),silent=TRUE)
  if(verbose){cat("the change in the log-likelihood is", Lout$value,"\n")}
  if(inherits(Lout,"try-error") || Lout$value > 199 ||
    Lout$value < -790) {
    if(verbose){
      cat("MLE could not be found. Degenerate!\n")
      cat("Nelder-Mead Log-likelihood ratio is ", Lout$value,"\n")
    }
    return(c(Inf, guess))
  }
  theta <- Lout$par
  names(theta) <- names(theta0)
# c0  <- llik.fun(theta=Lout$par, xobs=xobs,
#                 xsim=xsim, probs=probs,
#                 xsim.miss=xsim.miss, probs.miss=probs.miss,
#                 penalty=0.5, eta0=eta0, etamap=etamap)
  loglikelihood <- Lout$value
  names(loglikelihood) <- "delta.log.lik"

# loglikelihood
#
# Returns the change in the log-likelihood ratio if the dyad is toggled
# and the change in the coefficient
# list(coef=theta, 
#      loglikelihood=loglikelihood)
  c(loglikelihood, theta)
}
