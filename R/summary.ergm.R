summary.ergm <- function (object, ..., correlation=FALSE, covariance=FALSE)
{
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
  
  cat("\n==========================\n")
  cat("Summary of model fit\n")
  cat("==========================\n\n")
  
  cat("Formula:   ")
  print(object$formula)
  cat("\n")
  
#  cat("Estimate of network statistics for g:\n\n")
  {
#    print(object)

    digits = max(3, getOption("digits") - 3)
    cat ("Newton-Raphson iterations: ", object$iterations[1], "\n")
#
#   MSH changed to make clearer
#   original <- format(object$MCMCtheta, digits = digits)
    original <- format(object$theta.original, digits = digits)
    if(is.na(object$samplesize) & !all(object$theta1$independent)){
      cat ("\nPseudolikelihood Results:\n")
    }else{
      cat ("MCMC sample of size", object$samplesize, "\n")
      cat ("\nMonte Carlo MLE Results:\n")
    }
    
    if(!is.null(object$re)){ 
     if(!is.matrix(object$re)){
       cat ("\n Activity random effects:\n  Variances:\n")
       print(object$re)
     }else{
      cat ("\nSender and Receiver random effects:\n  Covariances:\n")
      print(object$re)       
      corr <- object$re[1,2]/sqrt(object$re[1,1]*object$re[2,2])
      corr <- max(min(1,corr),-1)
      cat (paste("\n  Correlation between sender and receiver:  ",
          round(corr,5)),"\n\n")
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
#   rdf <- dyads - length(object$coef)
    tval <- object$coef / asyse
    pval <- 2 * pt(q=abs(tval), df=rdf, lower.tail=FALSE)
    values <- format(object$coef,digits=digits)
    names <- names(values)
    names(values) <- NULL
    casyse<-format(asyse, digits=digits)
    cpval<-format(pval, digits=digits)
    cmc.se <- format(object$mc.se,digits=digits)

    cmc.se[object$offset] <- NA
    cpval[object$offset]  <- NA
    casyse[object$offset] <- NA

    count <- 1
    templist <- NULL
    while (count <= length(names))
      {
       templist <- append(templist,c(values[count],
            casyse[count],cpval[count],cmc.se[count]))
       count <- count+1
      }

    tempmatrix <- matrix(templist, ncol=4,byrow=TRUE)
    tempmatrix[,c(1:2,4)] <- format(tempmatrix[,c(1:2,4)], digits=digits, 
                                    print.gap=2)
    tempmatrix[,3] <- format.pval(as.numeric(tempmatrix[,3]), digits = 3, eps=1e-4)
    colnames(tempmatrix) <- c("estimate","s.e.","p-value","MCMC s.e.")
    rownames(tempmatrix) <- names
    print(tempmatrix, quote=FALSE)
  } 
  
  cat("\n")
  if(any(is.na(object$mc.se)) & !all(object$theta1$independent)){
   if(is.na(object$samplesize) & !all(object$theta1$independent)){
    warning("\n  The standard errors are based on naive pseudolikelihood and are suspect.\n")
   }else{
    warning("\n  The standard errors are suspect due to possible poor convergence.\n")
   }
  }
# if(is.na(object$samplesize)){
#   cat("Log Pseudo-likelihood: ",   object$mle.lik)
#   cat("\n    Pseudo-deviance:    ",-2*object$mle.lik,"on", rdf,"residual df.\n\n")
# }else{
#   cat("Log likelihood: ",object$mle.lik)
#   cat("\n    Deviance:    ",-2*object$mle.lik,"on", rdf, "residual df.\n\n")
# }

  if(!is.matrix(object$sample) & !all(object$theta1$independent)){
    cat(" ")
    cat(apply(cbind(paste(format(c("     Null", 
            "Residual", ""), width = 8, flag = ""), "Pseudo-deviance:"), 
            format(c(object$null.deviance,
                     -2*object$mle.lik, 
                     object$null.deviance+2*object$mle.lik),
                digits = 5), " on",
            format(c(dyads, rdf, df),
                digits = 5)," degrees of freedom\n"), 
            1, paste, collapse = " "),"\n")
  }else{
    cat(" ")
    cat(apply(cbind(paste(format(c("    Null", 
            "Residual", ""), width = 8), "Deviance:"), 
            format(c(object$null.deviance,
                     -2*object$mle.lik, 
                     object$null.deviance+2*object$mle.lik),
                digits = 5), " on",
            format(c(dyads, rdf, df),
                digits = 5)," degrees of freedom\n"), 
            1, paste, collapse = " "),"\n")
  }
# if(is.null(object$aic)){
   object$aic <- -2*object$mle.lik + 2*df
# }
# if(is.null(object$bic)){
   object$bic <- -2*object$mle.lik + log(dyads)*df
# }
  cat(paste("AIC:", format(object$aic, digits = 5), "  ", 
            "BIC:", format(object$bic, digits = 5), "\n", sep=" "))
  
  if (covariance == TRUE)
    {
      cat("Asymptotic covariance matrix:\n")
      print(asycov)
    }
  
#  cat("\nAsymptotic standard error vector:\n")
#  print(asyse)

  if (correlation == TRUE)
    {
      cat("\nAsymptotic correlation matrix:\n")
#     asycor <- diag(1/asyse) %*% asycov %*% diag(1/asyse)
      asycor <- asycov / crossprod(asyse)
      dimnames(asycor) <- dimnames(asycov)
      print(asycor)
    }
}
