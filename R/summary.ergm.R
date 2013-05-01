#  File R/summary.ergm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
###############################################################################
# The <summary.ergm> function prints a 'summary of model fit' table and returns
# the components of this table and several others listed below
#
# --PARAMETERS--
#   object     : an ergm object
#
#
# --IGNORED PARAMETERS--
#   ...        : used for flexibility
#   digits     : significant digits for the coefficients;default=
#                max(3,getOption("digits")-3), but the hard-coded value is 5
#   correlation: whether the correlation matrix of the estimated parameters
#                should be printed (T or F); default=FALSE
#   covariance : whether the covariance matrix of the estimated parameters
#                should be printed (T or F); default=FALSE
#   eps        : the numerical tolerance inputted to the native R function
#                <format.pval>; the code using 'eps' is now commented out;
#                default=.0001
#
# --RETURNED--
#   ans: a "summary.ergm" object as a list containing the following
#      formula         : object$formula
#      randomeffects   : object$re
#      digits          : the 'digits' inputted to <summary.ergm> or the default
#                        value (despite the fact the digits will be 5)
#      correlation     : the 'correlation' passed to <summary.ergm>
#      degeneracy.value: object$degenarcy.value
#      offset          : object$offset
#      drop            : object$drop
#      covariance      : the 'covariance' passed to <summary.ergm>
#      pseudolikelihood: whether pseudoliklelihood was used (T or F)
#      independence    : whether ?? (T or F)
#      iterations      : object$iterations
#      samplesize      : NA if 'pseudolikelihood'=TRUE, object$samplesize otherwise
#      message         : a message regarding the validity of the standard error
#                        estimates
#      aic             : the AIC goodness of fit measure
#      bic             : the BIC goodness of fit measure
#      coefs           : the dataframe of parameter coefficients and their
#                        standard erros and p-values
#      asycov          : the asymptotic covariance matrix
#      asyse           : the asymptotic standard error matrix
#      senderreceivercorrelation: 'randomeffects' if this is a matrix;
#                        otherwise, the correlation between sender and receiver??
#
################################################################################

summary.ergm <- function (object, ..., 
                          digits = max(3, getOption("digits") - 3),
                          correlation=FALSE, covariance=FALSE,
                          total.variation=TRUE,
                          eps=0.0001)
{
  control <- object$control
  pseudolikelihood <- object$estimate=="MPLE"
  independence <- is.dyad.independent(object)
  
  if(any(is.na(object$coef)) & !is.null(object$mplefit)){
     object$coef[is.na(object$coef)] <-
     object$mplefit$coef[is.na(object$coef)]
  }

  if(is.null(object$hessian) && is.null(object$covar)){
    object$covar <- diag(NA, nrow=length(object$coef))
  }
  
  nodes<- network.size(object$network)
  mc.se<- if(is.null(object$mc.se)) rep(NA, length(object$coef)) else object$mc.se
  

  if(is.null(object$covar)){
    asycov <- try(robust.inverse(-object$hessian), silent=TRUE)
    if(inherits(asycov,"try-error")){
      asycov <- diag(1/diag(-object$hessian))
    }
  }else{
    asycov <- object$covar
  }
  colnames(asycov) <- rownames(asycov) <- names(object$coef)
  
  asyse <- diag(asycov)
  if(total.variation && any(!is.na(mc.se))){
   asyse[!is.na(mc.se)] <- asyse[!is.na(mc.se)]+mc.se[!is.na(mc.se)]^2
  }
  asyse[asyse<0|is.infinite(object$coef)|object$offset] <- NA
  asyse <- sqrt(asyse)
  if(any(is.na(asyse)&!object$offset) & !is.null(object$mplefit)){
   if(is.null(object$mplefit$covar)){
    if(!is.null(object$mplefit$covar)){
     mpleasycov <- try(robust.inverse(-object$mplefit$hessian), silent=TRUE)
     if(inherits(mpleasycov,"try-error")){
      mpleasycov <- diag(1/diag(-object$mplefit$hessian))
     }
     asyse[is.na(asyse)] <- sqrt(diag(mpleasycov))[is.na(asyse)]
    }
   }else{
    mpleasycov <- object$mplefit$covar
    asyse[is.na(asyse)] <- sqrt(diag(mpleasycov))[is.na(asyse)]
   }
  }
  asyse <- matrix(asyse, ncol=length(asyse))
  colnames(asyse) <- colnames(asycov)

  ans <- list(formula=object$formula,
              digits=digits, correlation=correlation,
              degeneracy.value = object$degeneracy.value,
              offset = object$offset,
              drop = if(is.null(object$drop)) rep(0,length(object$offset)) else object$drop,
              estimable = if(is.null(object$estimable)) rep(TRUE,length(object$offset)) else object$estimable,
              covariance=covariance,
              pseudolikelihood=pseudolikelihood,
              independence=independence,
              estimate=object$estimate,
              control=object$control)
  
  ans$samplesize <- switch(object$estimate,
                           EGMME = if(!is.null(control$EGMME.main.method)) switch(control$EGMME.main.method,
                             `Gradient-Descent`=control$SA.phase3n,
                             stop("Unknown estimation method. This is a bug.")),
                           MPLE = NA,
                           MLE = if(!is.null(control$main.method)) switch(control$main.method,
                             `Stochastic-Approximation`=,
                             MCMLE=control$MCMC.samplesize,
                             `Robbins-Monro`=control$RM.phase3n,
                             `Stepping`=control$Step.MCMC.samplesize,
                             stop("Unknown estimation method. This is a bug.")),
                           stop("Unknown estimate type. This is a bug.")
                           )
                              

  ans$iterations <- switch(object$estimate,
                           EGMME = if(!is.null(control$EGMME.main.method)) switch(control$EGMME.main.method,
                             `Gradient-Descent`=NA,
                             stop("Unknown estimation method. This is a bug.")),
                           MPLE = NA,
                           MLE = if(!is.null(control$main.method)) switch(control$main.method,
                             `Stochastic-Approximation`=NA,
                             MCMLE=control$MCMLE.maxit,
                             `Robbins-Monro`=NA,
                             `Stepping`=NA,
                             stop("Unknown estimation method. This is a bug.")),
                           stop("Unknown estimate type. This is a bug.")
                           )
  
  nodes<- network.size(object$network)
  dyads<- network.dyadcount(object$network,FALSE)-network.edgecount(NVL(get.miss.dyads(object$constrained, object$constrained.obs),network.initialize(1)))
  df <- length(object$coef)

  rdf <- dyads - df
  tval <- object$coef / asyse
  pval <- 2 * pt(q=abs(tval), df=rdf, lower.tail=FALSE)

# Convert to % error
  if(any(!is.na(mc.se))){
   mc.se[!is.na(mc.se)] <- round(100*mc.se[!is.na(mc.se)]^2/asyse[!is.na(mc.se)]^2)
  }

  count <- 1
  templist <- NULL
  while (count <= length(names(object$coef)))
    {
     templist <- append(templist,c(object$coef[count],
          asyse[count],mc.se[count],pval[count]))
     count <- count+1
    }

  tempmatrix <- matrix(templist, ncol=4,byrow=TRUE)
  colnames(tempmatrix) <- c("Estimate", "Std. Error", "MCMC %", "p-value")
  rownames(tempmatrix) <- names(object$coef)

  devtext <- "Deviance:"
  if (object$estimate!="MPLE" || !independence || object$reference != as.formula(~Bernoulli)) {
    if (pseudolikelihood) {
      devtext <- "Pseudo-deviance:"
      ans$message <- "\nWarning:  The standard errors are based on naive pseudolikelihood and are suspect.\n"
    } 
    else if(object$estimate == "MLE" && any(is.na(mc.se)) && 
                      (!independence || control$force.main) ) {
      ans$message <- "\nWarning:  The standard errors are suspect due to possible poor convergence.\n"
    }
  } else {
    ans$message <- "\nFor this model, the pseudolikelihood is the same as the likelihood.\n"
  }
  mle.lik<-try(logLik(object,...), silent=TRUE)
  null.lik<-try(logLikNull(object,...), silent=TRUE)

  ans$null.lik.0 <- is.na(null.lik)

  if(!inherits(mle.lik,"try-error")){
  
    ans$devtable <- c("",apply(cbind(paste(format(c("    Null", "Residual"), width = 8), devtext), 
                                     format(c(if(is.na(null.lik)) 0 else -2*null.lik, -2*mle.lik), digits = digits), " on",
                                     format(c(dyads, rdf), digits = digits)," degrees of freedom\n"), 
                               1, paste, collapse = " "),"\n")
    
    ans$aic <- AIC(mle.lik)
    ans$bic <- BIC(mle.lik)
  }else ans$objname<-deparse(substitute(object))
  
  ans$coefs <- as.data.frame(tempmatrix)
  ans$asycov <- asycov
  ans$asyse <- asyse
  class(ans) <- "summary.ergm"
  ans
}

