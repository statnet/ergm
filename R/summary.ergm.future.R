Summary.ergm.future <- function (object, ..., 
                          digits = max(3, getOption("digits") - 3),
                          correlation=FALSE, covariance=FALSE,
                          eps=0.0001)
{
# separates out summary and print fns: MSH
  if(any(is.na(object$coef)) & !is.null(object$mplefit)){
     object$coef[is.na(object$coef)] <-
     object$mplefit$coef[is.na(object$coef)]
  }
  if(is.null(object$hessian) && is.null(object$covar)){
   return()
  }
  if(is.null(object$covar)){
   asycov <- try(robust.inverse(-object$hessian), silent=TRUE)
   if(inherits(asycov,"try-error")){
    asycov <- diag(1/diag(-object$hessian))
   }
  }else{
   asycov <- object$covar
  }
  rownames(asycov) <- names(object$coef)
  colnames(asycov) <- rownames(asycov)
  
  asyse <- diag(asycov)
  asyse[asyse<0] <- NA
  asyse <- sqrt(asyse)
  if(any(is.na(asyse)) & !is.null(object$mplefit)){
   if(is.null(object$mplefit$covar)){
    mpleasycov <- try(robust.inverse(-object$mplefit$hessian), silent=TRUE)
    if(inherits(mpleasycov,"try-error")){
     mpleasycov <- diag(1/diag(-object$mplefit$hessian))
    }
   }else{
    mpleasycov <- object$mplefit$covar
   }
   asyse[is.na(asyse)] <- sqrt(diag(mpleasycov))[is.na(asyse)]
  }
  asyse <- matrix(asyse, ncol=length(asyse))
  colnames(asyse) <- colnames(asycov)
  
#   MSH changed to make clearer
#   original <- format(object$MCMCtheta, digits = digits)
#   original <- format(object$theta.original, digits = digits)

  ans <- list(formula=object$formula, randomeffects=object$re,
              digits=digits, correlation=correlation,
              covariance=covariance,
              iterations=object$iterations[1])

  if(is.na(object$samplesize) & !all(object$theta1$independent)){
    ans$samplesize <- NA
  }else{
    ans$samplesize <-  object$samplesize
  }
    
  if(!is.null(object$re)){ 
   if(!is.matrix(object$re)){
    ans$senderreceivercorrelation<-object.re
   }else{
    corr <- object$re[1,2]/sqrt(object$re[1,1]*object$re[2,2])
    corr <- max(min(1,corr),-1)
    ans$senderreceivercorrelation<-corr
   }
  }

  nodes<- network.size(object$network)
  dyads<- network.dyadcount(object$network)
  if(!is.null(object$Z.mkl)){
    p <- ncol(object$Z.mkl)
  }else{
    p <- 0
  }
  if(!is.null(object$cluster)){
    df <- length(object$coef) + object$ngroups*(p+2) - 1 # ng-1 + ng *p + ng
  }else{
    df <- length(object$coef) + (nodes - (p + 1)/2) * p
  }
  rdf <- dyads - df
  tval <- object$coef / asyse
  pval <- 2 * pt(q=abs(tval), df=rdf, lower.tail=FALSE)

  count <- 1
  templist <- NULL
  while (count <= length(names(object$coef)))
    {
     templist <- append(templist,c(object$coef[count],
          asyse[count],object$mc.se[count],pval[count]))
     count <- count+1
    }

  tempmatrix <- matrix(templist, ncol=4,byrow=TRUE)
  colnames(tempmatrix) <- c("Estimate", "Std. Error", "MCMC s.e.", "p-value")
  rownames(tempmatrix) <- names(object$coef)
  
if(any(is.na(object$mc.se)) & !all(object$theta1$independent)){
 if(is.na(object$samplesize) & !all(object$theta1$independent)){
  ans$warning <- "\n  The standard errors are based on naive pseudolikelihood and are suspect.\n"
 }else{
  ans$warning <- "\n  The standard errors are suspect due to possible poor convergence.\n"
 }
}

  if(!is.matrix(object$sample) & !all(object$theta1$independent)){
    ans$devtable <- c("",apply(cbind(paste(format(c("    Null", 
            "Residual", ""), width = 8, flag = ""), "Pseudo-deviance:"), 
            format(c(object$null.deviance,
                     -2*object$mle.lik, 
                     object$null.deviance+2*object$mle.lik),
                digits = 5), " on",
            format(c(dyads, rdf, df),
                digits = 5)," degrees of freedom\n"), 
            1, paste, collapse = " "),"\n")
  }else{
    ans$devtable <- c("",apply(cbind(paste(format(c("   Null", 
            "Residual", ""), width = 8), "Deviance:"), 
            format(c(object$null.deviance,
                     -2*object$mle.lik, 
                     object$null.deviance+2*object$mle.lik),
                digits = 5), " on",
            format(c(dyads, rdf, df),
                digits = 5)," degrees of freedom\n"), 
            1, paste, collapse = " "),"\n")
  }
  ans$aic <- -2*object$mle.lik + 2*df
  ans$bic <- -2*object$mle.lik + log(dyads)*df
  
  ans$coefs <- tempmatrix
  ans$asycov <- asycov
  ans$asyse <- asyse
  class(ans) <- "summary.ergm"
  ans
}

print.summary.ergm <- function (x, 
              digits = max(3, getOption("digits") - 3),
              correlation=FALSE, covariance=FALSE,
              signif.stars= getOption("show.signif.stars"),
              eps=0.0001, ...)
{
  cat("\n==========================\n")
  cat("Summary of model fit\n")
  cat("==========================\n\n")
  
  cat("Formula:   ")
  print(x$formula)
  cat("\n")
  
  cat ("Newton-Raphson iterations: ", x$iterations, "\n")
  if(is.na(x$samplesize)){
    cat ("\nPseudolikelihood Results:\n")
  }else{
    cat ("MCMC sample of size", x$samplesize, "\n")
    cat ("\nMonte Carlo MLE Results:\n")
  }
    
  if(!is.null(x$randomeffects)){ 
     if(!is.matrix(x$randomeffects)){
       cat ("\n Activity random effects:\n  Variances:\n")
       print(x$randomeffects)
     }else{
      cat ("\nSender and Receiver random effects:\n  Covariances:\n")
      print(x$randomeffects)       
      corr <- x$randomeffects[1,2]/sqrt(x$randomeffects[1,1]*x$randomeffects[2,2])
      corr <- max(min(1,corr),-1)
      cat (paste("\n  Correlation between sender and receiver:  ",
          round(corr,5)),"\n\n")
     }
  }

  printCoefmat(x$coefs, digits=digits, signif.stars=signif.stars,
               P.values=TRUE, has.Pvalue=TRUE, na.print="NA",
               eps=eps, ...)
  
  if(!is.null(x$warning)){ warning(x$warning) }

  cat("\n")
  cat(x$devtable)

  cat(paste("AIC:", format(x$aic, digits = 5), "  ", 
            "BIC:", format(x$bic, digits = 5), "\n", sep=" "))
  

  if (covariance == TRUE)
    {
      cat("Asymptotic covariance matrix:\n")
      print(x$asycov)
    }
  
  if (correlation == TRUE)
    {
      cat("\nAsymptotic correlation matrix:\n")
      asycor <- x$asycov / crossprod(x$asyse)
      dimnames(asycor) <- dimnames(x$asycov)
      print(asycor)
    }
  
  invisible(x)
}
