###############################################################################
# The <print.summary.ergm> function prints a subset of the information given
# by <summary.ergm>
#
# --PARAMETERS--
#   x           : a "summary.ergm" object, as returned by <summary.ergm>
#   digits      : the number of significant digits for the coefficients;
#                 default=max(3, getOption("digits")-3)
#   correlation : whether the correlation matrix of the estimated parameters
#                 should be printed (T or F); default=FALSE
#   covariance  : whether the covariance matrix of the estimated parameters
#                 should be printed (T or F); default=FALSE
#   signif.stars: whether stars are to be printed on summary tables of
#                 coefficients (T or F); default=getOption("show.signif.stars")
#   eps.Pvalue  : the tolerance to be passed to the R <printCoefmat> function;
#                 default=.0001
#   ...         : additional parameters to be passed to <printCoefmat> 
#
# --RETURNED--
#   x
#
###############################################################################

print.summary.ergm <- function (x, 
              digits = max(3, getOption("digits") - 3),
              correlation=FALSE, covariance=FALSE,
              signif.stars= getOption("show.signif.stars"),
              eps.Pvalue=0.0001, print.header=TRUE, print.formula=TRUE, print.fitinfo=TRUE, print.coefmat=TRUE, print.message=TRUE, print.deviances=TRUE, print.drop=TRUE, print.offset=TRUE, print.degeneracy=TRUE,...){
  if(print.header){
    cat("\n==========================\n")
    cat("Summary of model fit\n")
    cat("==========================\n\n")
  }

  if(print.formula){
    cat("Formula:   ")
    print(x$formula)
    cat("\n")
  }

  if(print.fitinfo){
    if (!is.null(x$iterations)) {
      cat("Newton-Raphson iterations: ", x$iterations, "\n")
    }
  
    if(x$pseudolikelihood){
      if (x$independence) {
        cat ("\nMaximum Likelihood Results:\n")
      } else {
        cat ("\nMaximum Pseudolikelihood Results:\n")
      }
    }else{
      cat ("MCMC sample of size", x$samplesize, "\n")
      cat ("\nMonte Carlo MLE Results:\n")
    }
  }

  if(print.coefmat){
    printCoefmat(x$coefs, digits=digits, signif.stars=signif.stars,
                 P.values=TRUE, has.Pvalue=TRUE, na.print="NA",
                 eps.Pvalue=eps.Pvalue, ...)
  }

  if(print.message){
    if(!is.null(x$message)){ 
      cat(x$message)
    }
    cat("\n")
  }

  if(print.deviances){
    if(!is.null(x$devtable)){
      cat(x$devtable)
      
      cat(paste("AIC:", format(x$aic, digits = 5), "  ", 
                "BIC:", format(x$bic, digits = 5), "\n", sep=" "))
    } else cat(nologLik.message(x$objname))
  }

  if(print.drop){
    if(any(x$drop!=0)){
      cat("\n Warning: The following terms have infinite coefficient estimates:\n  ")
      cat(rownames(x$coefs)[x$drop!=0], "\n\n")
    }
  }

  if(print.offset){
    if(any(x$offset&!x$drop!=0)){
      cat("\n Warning: The following terms are fixed by offset and are not estimated:\n  ")
      cat(rownames(x$coefs)[x$offset & !x$drop!=0], "\n\n")
    }
  }

  if(print.degeneracy){
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
  }

  if(covariance == TRUE){
    cat("Asymptotic covariance matrix:\n")
    print(x$asycov)
  }
  
  if(correlation == TRUE){
    cat("\nAsymptotic correlation matrix:\n")
    asycor <- x$asycov / crossprod(x$asyse)
    dimnames(asycor) <- dimnames(x$asycov)
    print(asycor)
  }
  
  invisible(x)
}
