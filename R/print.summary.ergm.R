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
  if(x$pseudolikelihood){
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
  
  if(!is.null(x$warning)){ 
     cat("\n Warning:\n")
     cat(x$warning)
  }

  cat("\n")
  cat(x$devtable)

  cat(paste("AIC:", format(x$aic, digits = 5), "  ", 
            "BIC:", format(x$bic, digits = 5), "\n", sep=" "))
  

  if(any(x$drop)){
    cat("\n Warning:\n")
    for(i in names(x$coefs[x$offset,1])){
     cat(paste("  The term",i,
     "is degenerate and has an infinite coefficient estimate.\n",
      sep=" "))
    }
  }

  if(any(x$offset&!x$drop)){
    cat("\n Warning:\n")
    for(i in names(x$coefs[x$offset,1])){
    cat(paste("  The term",i,
     "has been offset and was not estimated from the data.\n",
      sep=" "))
    }
  }

  if(!is.null(x$degeneracy.value) && !is.na(x$degeneracy.value)){
   if(is.infinite(x$degeneracy.value)){
    cat("\n Warning: The diagnostics indicate that the model is very unstable.\n   They suggest that the model is near degenerate,\n   and that the numerical summaries are suspect.\n")
   }else{
    if(x$degeneracy.value > 1){
      cat("The instability of the model is: ",
        format(x$degeneracy.value, digits=2),"\n")
      cat("Instabilities greater than 1 suggest the model is near degenerate.\n")
    }
   }
  }

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
