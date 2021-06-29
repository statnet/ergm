#  File R/print.summary.ergm.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################

#' @rdname summary.ergm
#' @order 2
#' 
#' @param x object of class `summary.ergm` returned by [summary.ergm()].
#' @param digits significant digits for coefficients
#' @param signif.stars whether to print dots and stars to signify
#'   statistical significance. See [print.summary.lm()].
#' @param eps.Pvalue \eqn{p}-values below this level will be printed
#'   as "<`eps.Pvalue`".
#' @param
#'   print.formula,print.fitinfo,print.coefmat,print.message,print.deviances,print.drop,print.offset,print.call
#'   which components of the fit summary to print.
#'   
#' @details The default printout of the summary object contains the
#'   call, number of iterations used, null and residual deviances, and
#'   the values of AIC and BIC (and their MCMC standard errors, if
#'   applicable). The coefficient table contains the following
#'   columns:
#'   
#'   - `Estimate`, `Std. Error` - parameter estimates and their standard errors
#'   - `MCMC %` - if `total.variation=TRUE` (default) the percentage of standard
#'   error attributable to MCMC estimation process rounded to an integer. See
#'   also [vcov.ergm()] and its `sources` argument.
#'   - `z value`, `Pr(>|z|)` - z-test and p-values
#'    
#' @export
print.summary.ergm <- function (x, 
              digits = max(3, getOption("digits") - 3),
              correlation=x$correlation, covariance=x$covariance,
              signif.stars= getOption("show.signif.stars"),
              eps.Pvalue=0.0001, print.formula=FALSE, print.fitinfo=TRUE, print.coefmat=TRUE, print.message=TRUE, print.deviances=TRUE, print.drop=TRUE, print.offset=TRUE, print.call=TRUE,...){
  
  control <- x$control

  # The following code is based on stats:::print.lm(), but there really isn't another concise way to do this:
  if(print.call && !is.null(x$call)) cat("Call:\n", paste(deparse(x$call), sep="\n", collapse="\n"), "\n\n", sep="")

  if(print.formula) cat("Formula:\n", paste(deparse(x$formula), sep="\n", collapse="\n"), "\n\n", sep="")

  if(print.fitinfo){
    ## if (!is.null(x$iterations)) {
    ##   cat("Iterations: ", x$iterations, "\n")
    ## }

    cat(paste0(x$estimate.desc, " Results:\n\n"))
  }

  if(print.coefmat){
    printCoefmat(coef(x), digits=digits, signif.stars=signif.stars,
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
    if(is.null(x$devtable)) message(NO_LOGLIK_MESSAGE)
    else if(length(x$devtable)>1 || !is.na(x$devtable)){
      cat(c("",apply(cbind(paste(format(c("    Null", "Residual"), width = 8), x$devtext), 
                                     format(x$devtable[,1], digits = digits), " on",
                                     format(x$devtable[,2], digits = digits)," degrees of freedom\n"), 
                               1, paste, collapse = " "),"\n"))

      if(x$null.lik.0) writeLines(c(strwrap(paste("Note that the null model likelihood and deviance are defined to be 0.", NO_NULL_IMPLICATION)),''))
      
      cat(paste0("AIC: ", format(x$aic, digits = digits), "  ",
                 "BIC: ", format(x$bic, digits = digits), "  ",
                 "(Smaller is better. MC Std. Err. = ", format(sqrt(NVL(attr(x$aic,"vcov"),0)), digits=digits), ")", "\n"))
    }
  }

  if(print.drop){
    if(any(x$drop!=0)){
      cat("\n Warning: The following terms have infinite coefficient estimates:\n  ")
      cat(rownames(coef(x))[x$drop!=0], "\n")
    }
    if(any(!x$estimable)){
      cat("\n Warning: The following terms could not be estimated because they conflicted with the sample space constraint:\n  ")
      cat(rownames(coef(x))[!x$estimable], "\n")
    }
  }

  if(print.offset){
    if(any(x$offset & x$drop==0 & x$estimable)){
      cat("\n The following terms are fixed by offset and are not estimated:\n  ")
      cat(rownames(coef(x))[x$offset & x$drop==0 & x$estimable], "\n\n")
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
