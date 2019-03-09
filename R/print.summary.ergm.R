#  File R/print.summary.ergm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################
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

#' @rdname summary.ergm
#' @param x object of class `summary.ergm` returned by
#'   [summary.ergm()].
#' @param signif.stars whether to print dots and stars to signify
#'   statistical significance. See [print.summary.lm()].
#' @param eps.Pvalue \eqn{p}-values below this level will be printed
#'   as "<`eps.Pvalue`".
#' @param
#'   print.header,print.formula,print.fitinfo,print.coefmat,print.message,print.deviances,print.drop,print.offset,print.degeneracy
#'   which components of the fit summary to print.
#' @export
print.summary.ergm <- function (x, 
              digits = max(3, getOption("digits") - 3),
              correlation=FALSE, covariance=FALSE,
              signif.stars= getOption("show.signif.stars"),
              eps.Pvalue=0.0001, print.header=TRUE, print.formula=TRUE, print.fitinfo=TRUE, print.coefmat=TRUE, print.message=TRUE, print.deviances=TRUE, print.drop=TRUE, print.offset=TRUE, print.degeneracy=TRUE,...){
  
  control <- x$control
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
      cat("Iterations: ", x$iterations, "\n")
    }

    switch(x$estimate,
           MPLE = if (x$independence) {
             cat("\nMaximum Likelihood Results:\n")
           } else {
             cat("\nMaximum Pseudolikelihood Results:\n")
           },
           CD = cat("\nContrastive Divergence results:\n"),
           MLE = NVL3(control$main.method, switch(.,
             MCMLE = cat("\nMonte Carlo MLE Results:\n"),
             `Stochastic-Approximation`=cat("\nMonte Carlo MLE Results:\n"),
             `Robbins-Monro`=cat("\nRobbins-Monro MLE Results:\n"),
             `Stepping`=cat("\n Stepping MLE Results:\n"),
             stop("Unknown estimation method. This is a bug."))),
           EGMME = NVL3(control$EGMME.main.method, switch(.,
             `Gradient-Descent`=cat("\nEquilibrium Generalized Method of Moments Results:\n"),
             stop("Unknown estimation method. This is a bug."))),
           stop("Unknown estimate type. This is a bug.")
           )
  }

  if(print.coefmat){
    printCoefmat(x$coefficients, digits=digits, signif.stars=signif.stars,
                 P.values=TRUE, has.Pvalue=TRUE, na.print="NA",
                 eps.Pvalue=eps.Pvalue, cs.ind=1:2, tst.ind=4L,...)
  }

  if(print.message){
    if(!is.null(x$message)){ 
      cat(x$message)
    }
    cat("\n")
  }

  if(print.deviances){
    if(!is.null(x$devtable)){
      cat(c("",apply(cbind(paste(format(c("    Null", "Residual"), width = 8), x$devtext), 
                                     format(x$devtable[,1], digits = digits), " on",
                                     format(x$devtable[,2], digits = digits)," degrees of freedom\n"), 
                               1, paste, collapse = " "),"\n"))

      if(x$null.lik.0) writeLines(c(strwrap(paste("Note that the null model likelihood and deviance are defined to be 0.", NO_NULL_IMPLICATION)),''))
      
      cat(paste("AIC:", format(x$aic, digits = digits), "  ", 
                "BIC:", format(x$bic, digits = digits), "  ",
                "(Smaller is better.)", "\n", sep=" "))
    } else cat(nologLik.message(x$objname))
  }

  if(print.drop){
    if(any(x$drop!=0)){
      cat("\n Warning: The following terms have infinite coefficient estimates:\n  ")
      cat(rownames(x$coefficients)[x$drop!=0], "\n")
    }
    if(any(!x$estimable)){
      cat("\n Warning: The following terms could not be estimated because they conflicted with the sample space constraint:\n  ")
      cat(rownames(x$coefficients)[!x$estimable], "\n")
    }
  }

  if(print.offset){
    if(any(x$offset & x$drop==0 & x$estimable)){
      cat("\n The following terms are fixed by offset and are not estimated:\n  ")
      cat(rownames(x$coefficients)[x$offset & x$drop==0 & x$estimable], "\n\n")
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
