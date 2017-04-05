#  File ergm/R/summary.ergm.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
###############################################################################
# The <summary.ergm> function prints a 'summary of model fit' table and returns
# the components of this table and several others listed below
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
  dyads<- network.dyadcount(object$network)
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
                             `Stochastic-Approximation`=control$SA.phase3n,
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
                             `Stochastic-Approximation`=NA,
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
  dyads<- network.dyadcount(object$network)
  df <- length(object$coef)

  rdf <- dyads - df
  tval <- object$coef / asyse
  pval <- 2 * pt(q=abs(tval), df=rdf, lower.tail=FALSE)

# Convert to % error
  if(any(!is.na(mc.se))){
   mc.se[!is.na(mc.se)] <- round(100*mc.se[!is.na(mc.se)]*mc.se[!is.na(mc.se)]/(asyse[!is.na(mc.se)]*asyse[!is.na(mc.se)]))
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
  if (object$estimate!="MPLE" || !independence) {
    if (pseudolikelihood) {
      devtext <- "Pseudo-deviance:"
      ans$message <- "\nWarning:  The standard errors are based on naive pseudolikelihood and are suspect.\n"
    } 
    else if(object$estimate == "MLE" && any(is.na(mc.se)) && 
                      (!independence || control$force.main) ) {
      ans$message <- "\nWarning:  The standard errors are suspect due to possible poor convergence.\n"
    }else if(object$estimate == "EGMME"){
      ans$message <- "\nWarning:  The standard errors do not incorporate uncertainty due to the \"noisy\" estimation procedure.\n"
    }
  } else {
    ans$message <- "\nFor this model, the pseudolikelihood is the same as the likelihood.\n"
  }
  llk<-try(logLik(object,...), silent=TRUE)

  if(!inherits(llk,"try-error")){
  
    ans$devtable <- c("",apply(cbind(paste(format(c("   Null", 
                                                    "Residual", ""), width = 8), devtext), 
                                     format(c(object$null.deviance,
                                              -2*llk, 
                                              object$null.deviance+2*llk),
                                            digits = 5), " on",
                                     format(c(dyads, rdf, df),
                                            digits = 5)," degrees of freedom\n"), 
                               1, paste, collapse = " "),"\n")
    
    
    ans$aic <- AIC(llk)
    ans$bic <- BIC(llk)
  }else ans$objname<-deparse(substitute(object))
  
  ans$coefs <- as.data.frame(tempmatrix)
  ans$asycov <- asycov
  ans$asyse <- asyse
  class(ans) <- "summary.ergm"
  ans
}

