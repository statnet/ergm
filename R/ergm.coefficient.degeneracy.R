#  File R/ergm.coefficient.degeneracy.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2015 Statnet Commons
#######################################################################
#######################################################################
# The <ergm.coefficient.degeneracy> function
#
# --PARAMETERS--
#   object   : an ergm object
#   control  : a control list, whose only recognized component is:
#     trustregion:  loglikelihoods above ('trustregion' - 0.1) obtain
#                    -Inf degeneracy.values
#              default= the list returned by <control.ergm>, which 
#              yields 'trustregion'=20 
#   test.only: whether the instability of the model should not be
#              printed (T or F); if FALSE, the instability is printed;
#              if TRUE, the instability is not printed, unless the
#              returned 'degeneracy.type' is > 1; default=FALSE
#   verbose  : whether the returned 'degeneracy.type' should be printed
#              (T or F); default=FALSE
#
# --IGNORED--
#   fast     : ?? ; default=TRUE 
#
# --RETURNED--
#   the orginal ergm object with 2 new additional components:
#      degeneracy.value: the degeneracy value, which may take on one of
#                        three values:
#            NA - if the ergm object's glm mplefit doesn't exist and 
#                 cannot be fitted
#            Inf- if the ergm's log likelihood is greater than
#                'trustregion' - 0.1   or  if the non-MCMC fit of the
#                 ergm doesn't have "glm" in its class
#            the max 'degeneracy.type' value - in all other cases
#      degeneracy.type :
#            NULL - if the ergm object's glm mplefit doesn't exist and
#                   cannot be fitted
#            the loglikelihood of the MCMC fit - if an MCMC fit
# 
#######################################################################

ergm.coefficient.degeneracy <- function(object, 
                          control=object$control,
                          fast=TRUE,
                          test.only=FALSE,
                          verbose=FALSE) {
  check.control.class(control, "ergm")
  if(!is.ergm(object)){
    stop("A ergm object argument must be given.")
  }
  if(is.matrix(object$sample)){
   if(is.null(object$mplefit$glm)){
    current.warn <- options()$warn
    options(warn=-1)
    fit <- try(ergm(object$formula, estimate="MPLE"), silent=TRUE)
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
   if(object$loglikelihood>control$MCMLE.trustregion-0.1){
    object$degeneracy.value <- Inf
   }else{
    changeobs <- object$MCMCtheta-object$MCMCtheta
    wgts <- object$mplefit$glm$prior.weights
    object$degeneracy.type <- rep(Inf, length(object$MCMCtheta))
    names(object$degeneracy.type) <- names(object$MCMCtheta)
    for(i in seq(along=object$MCMCtheta)){
     etamapi <- object$etamap
     etamapi$offsettheta[-i] <- TRUE
     degeneracy.value <- try(
      ergm.compute.degeneracy(changeobs,
      init=object$MCMCtheta, etamap=etamapi, 
      statsmatrix=object$sample[,!object$etamap$offsettheta,drop=FALSE],
      trustregion=control$MCMLE.trustregion, guess=object$MCMCtheta[i]),silent=TRUE)
     if(!inherits(degeneracy.value,"try-error")){
      object$degeneracy.type[i] <- degeneracy.value[1]
     }
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
    object$degeneracy.type <- object$degeneracy.type[,1]
    object$degeneracy.value <- max(object$degeneracy.type,na.rm=TRUE)
   }else{
    object$degeneracy.value <- Inf
    object$degeneracy.type <- NULL
   }
  }
  if(any(object$degeneracy.type>control$MCMLE.trustregion-0.1)){
   object$degeneracy.type[object$degeneracy.type>control$MCMLE.trustregion-0.1] <- Inf
  }
  if(!is.null(object$degeneracy.type)){
   object$degeneracy.value <- max(object$degeneracy.type,na.rm=TRUE)
  }
  if(any(is.infinite(object$degeneracy.type))){
   cat("\n Warning: The diagnostics indicate that the model is very unstable.\n   They suggest that the model is degenerate,\n   and that the numerical summaries are suspect.\n")
   cat("\n The suspect coefficients are:\n      ")
   cat(paste(names(object$degeneracy.type)[is.infinite(object$degeneracy.type)],"\n     "))
  }else{
    if(!test.only || any(object$degeneracy.type > 1)){
     cat("The instability of the model is: ",
        format(object$degeneracy.value, digits=2),"\n")
    }
    if(any(object$degeneracy.type > 1)){
      cat("Instabilities greater than 1 suggest the model is degenerate.\n")
    }
  }
  if(verbose){
    print(object$degeneracy.type)
  }
  return(invisible(object))
}
