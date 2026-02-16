#  File R/summary.ergm.R in package ergm, part of the Statnet suite of packages
#  for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################

#' Summarizing ERGM Model Fits
#'
#' [base::summary()] method for [ergm()] fits.
#'
#' @param object an object of class `ergm`, usually, a result of a call to
#'   [ergm()].
#' @param correlation logical; if `TRUE`, the correlation matrix of the
#'   estimated parameters is returned and printed.
#' @param covariance logical; if `TRUE`, the covariance matrix of the estimated
#'   parameters is returned and printed.
#' @param total.variation logical; if `TRUE`, the standard errors reported in
#'   the `Std. Error` column are based on the sum of the likelihood variation
#'   and the MCMC variation. If `FALSE` only the likelihood variation is used.
#'   The \eqn{p}-values are based on this source of variation.
#' @param ... For [summary.ergm()] additional arguments are passed to
#'   [logLik.ergm()]. For [print.summary.ergm()], to [stats::printCoefmat()].
#'
#' @details [summary.ergm()] tries to be smart about formatting the
#' coefficients, standard errors, etc.
#'
#' @return The returned object is a list of class "ergm.summary" with
#'   the following elements:
#'   
#' \item{formula}{ERGM model formula}
#' \item{call}{R call used to fit the model}
#' \item{correlation, covariance}{whether to print correlation/covariance matrices of the estimated parameters}
#' \item{pseudolikelihood}{was the model estimated with MPLE}
#' \item{independence}{is the model dyad-independent}
#' \item{control}{the [control.ergm()] object used}
#' \item{samplesize}{MCMC sample size}
#' \item{message}{optional message on the validity of the standard error estimates}
#' \item{null.lik.0}{It is `TRUE` of the null model likelihood has not been calculated. See [logLikNull()]}
#' \item{devtext, devtable}{Deviance type and table}
#' \item{aic, bic}{values of AIC and BIC}
#' \item{coefficients}{matrices with model parameters and associated statistics}
#' \item{asycov}{asymptotic covariance matrix}
#' \item{asyse}{asymptotic standard error matrix}
#' \item{offset, drop, estimate, iterations, mle.lik, null.lik}{
#' see documentation of the object returned by [ergm()]
#' }
#'
#' @seealso The model fitting function [ergm()], [print.ergm()], and
#'   [base::summary()]. Function [stats::coef()] will extract the matrix of
#'   coefficients with standard errors, t-statistics and p-values.
#'
#' @keywords regression models
#' @examples
#'
#'  data(florentine)
#'
#'  x <- ergm(flomarriage ~ density)
#'  summary(x)
#'
#' @export
summary.ergm <- function (object, ..., 
                          correlation=FALSE, covariance=FALSE,
                          total.variation=TRUE)
{
  # Warn if the object was produced by an earlier version of ergm.
  myver <- packageVersion("ergm")
  objver <- NVL(object$ergm_version, as.package_version("3.9.4")) # 3.9.4 was the last version that didn't store the version information.
  nextver <- as.package_version(paste(objver$major, objver$minor+1, sep="."))
  if(objver < paste(myver$major, myver$minor, sep=".")){
    warn(paste0("This object was fit with ", sQuote("ergm"), " version ", objver, " or earlier. Summarizing it with version ", nextver, " or later may return incorrect results or fail."))
  }

  if("digits" %in% ...names()) warn("summary.ergm() no lnger takes a digits= argument.")
  control <- object$control
  pseudolikelihood <- object$estimate=="MPLE"
  independence <- NVL(object$MPLE_is_MLE, is.dyad.independent(object))
  coef <- coef(object)

  ans <- list(formula=object$formula,
              call=object$call,
              correlation=correlation,
              offset = object$offset,
              drop = NVL(object$drop, rep(0,length(object$offset))),
              estimable = NVL(object$estimable, rep(TRUE,length(object$offset))),
              covariance=covariance,
              pseudolikelihood=pseudolikelihood,
              independence=independence,
              estimate=object$estimate,
              estimate.desc=object$estimate.desc,
              control=object$control,
              lindep = object$lindep)
  
  asycov <- vcov(object, sources=if(total.variation) "all" else "model")
  asyse <- sqrt(diag(asycov))
  # Convert to % error  
  est.se <- sqrt(diag(vcov(object, sources="estimation")))
  mod.se <- sqrt(diag(vcov(object, sources="model")))
  tot.se <- sqrt(diag(vcov(object, sources="all")))
  est_pct <- ifelse(est.se > 0, round(100 * (tot.se - mod.se) / tot.se), 0)

  zval <- coef / asyse
  pval <- 2 * pnorm(q=abs(zval), lower.tail=FALSE)
  
  count <- 1
  coefmat <- cbind(
    `Estimate` = coef,
    `Std. Error` = asyse,
    `MCMC %` = est_pct,
    `z value` = zval,
    `Pr(>|z|)` = pval)

  rownames(coefmat) <- param_names(object)

  devtext <- "Deviance:"
  if (object$estimate!="MPLE" || !independence || object$reference != as.formula(~Bernoulli)) {
    if (pseudolikelihood && control$MPLE.covariance.method=="invHess") {
      devtext <- "Pseudo-deviance:"
      ans$message <- "\nWarning:  The standard errors are based on naive pseudolikelihood and are suspect. Set control.ergm$MPLE.covariance.method='Godambe' for a simulation-based approximation of the standard errors.\n"
    } 
    else if(object$estimate == "MLE" && any(is.na(est.se) & !ans$offset & !ans$drop==0 & !ans$estimable) && 
                      (!independence || control$force.main) ) {
      ans$message <- "\nWarning:  The standard errors are suspect due to possible poor convergence.\n"
    }
  } else {
    ans$message <- "\nFor this model, the pseudolikelihood is the same as the likelihood.\n"
  }
  mle.lik<-try(logLik(object,...), silent=TRUE)

  if(inherits(mle.lik,"try-error")) ans$objname<-deparse(substitute(object))
  else if(!is.na(mle.lik)){
    # Only evaluate the null likelihood if the MLE likelihood is defined.
    null.lik<-try(logLikNull(object,...), silent=TRUE)
    ans$null.lik.0 <- is.na(null.lik)

    df <- length(coef)
    dyads <- nobs(object)
    rdf <- dyads - df
    ans$devtable <- matrix(c(if(is.na(null.lik)) 0 else -2*null.lik, -2*mle.lik,
                             c(dyads, rdf)), 2,2, dimnames=list(c("Null","Residual"),
                                                                c("Resid. Dev", "Resid. Df")))
    ans$devtext <- devtext

    ans$aic <- AIC(object)
    ans$bic <- BIC(object)
    ans$mle.lik <- ERRVL(mle.lik, NA)
    ans$null.lik <- ERRVL(null.lik, NA)
  }else ans$devtable <- NA

  ans$coefficients <- coefmat
  ans$asycov <- asycov
  ans$asyse <- asyse
  class(ans) <- "summary.ergm"
  ans
}

#' @rdname summary.ergm
#' 
#' @param x object of class `summary.ergm` returned by [summary.ergm()].
#' @param digits significant digits for coefficients
#' @param signif.stars whether to print dots and stars to signify
#'   statistical significance. See [print.summary.lm()].
#' @param eps.Pvalue \eqn{p}-values below this level will be printed
#'   as "<`eps.Pvalue`".
#' @param
#'   print.formula,print.fitinfo,print.coefmat,print.message,print.deviances,print.drop,print.lindep,print.offset,print.call
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
              eps.Pvalue = 0.0001, print.formula = FALSE, print.fitinfo = TRUE,
              print.coefmat = TRUE, print.message = TRUE, print.deviances = TRUE,
              print.drop = TRUE, print.offset = TRUE, print.call = TRUE,
              print.lindep = TRUE, ...) {
  
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

  short_list <- function(x) paste(x, collapse = ", ")
  wrap_list <- function(x, indent = 3, exdent = 3) strwrap(x, indent = indent, exdent = exdent)

  if(print.offset){
    if(any(x$offset & x$drop==0 & x$estimable)){
      cat("\n The following terms are fixed by offset and are not estimated:\n")
      rownames(coef(x))[x$offset & x$drop==0 & x$estimable] |>
        quote_var_name() |>
        short_list() |>
        wrap_list() |>
        paste("\n") |>
        cat("\n")
    }
  }

  warnings <- c()

  wrap_head <- function(x) strwrap(x, exdent = 3, initial = " * ")

  if(print.drop){
    if (any(x$drop != 0))
      warnings <- c(warnings, "",
                    wrap_head("The following terms have infinite coefficient estimates due to an extreme sufficient statistic:"),
                    "",
                    rownames(coef(x))[x$drop != 0] |>
                    quote_var_name() |>
                    short_list() |>
                    wrap_list())

    if (any(!x$estimable))
      warnings <- c(warnings, "",
                    wrap_head("The following terms could not be estimated because they conflicted with the sample space constraint:"),
                    "",
                    rownames(coef(x))[!x$estimable] |>
                    quote_var_name() |>
                    short_list() |>
                    wrap_list())
  }

  if (print.lindep) {
    if(nrow(x$lindep) %||% 0)
      warnings <- c(warnings, "",
                    wrap_head("The following linear dependence has been detected among the model statistics and/or estimating functions:"),
                    "",
                    sapply(format(x$lindep), function(z) c(wrap_list(z, exdent = 5), "")) |> unlist()
                    )
  }

  if(length(warnings)) {
    cat("Warnings:\n")
    cat(paste0(warnings, "\n", collapse = ""))
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
