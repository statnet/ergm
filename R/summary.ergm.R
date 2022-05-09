#  File R/summary.ergm.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################



#' Summarizing ERGM Model Fits
#'
#' [base::summary()] method for [ergm()] fits.
#'
#' @order 1
#'
#' @param object an object of class "ergm", usually, a result of a call to
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

  if("digits" %in% names(list(...))) warn("summary.ergm() no lnger takes a digits= argument.")
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
              control=object$control)
  
  asycov <- vcov(object, sources=if(total.variation) "all" else "model")
  asyse <- sqrt(diag(asycov))
  # Convert to % error  
  est.se <- sqrt(diag(vcov(object, sources="estimation")))
  mod.se <- sqrt(diag(vcov(object, sources="model")))
  tot.se <- sqrt(diag(vcov(object, sources="all")))
  est.pct <- rep(NA,length(est.se))
  if(any(!is.na(est.se))){
    # We want (sqrt(V.model + V.MCMC)-sqrt(V.model))/sqrt(V.model + V.MCMC) * 100%,
    est.pct[!is.na(est.se)] <- ifelse(est.se[!is.na(est.se)]>0, round(100*(tot.se[!is.na(est.se)]-mod.se[!is.na(est.se)])/tot.se[!is.na(est.se)]), 0)
  }

  zval <- coef / asyse
  pval <- 2 * pnorm(q=abs(zval), lower.tail=FALSE)
  
  count <- 1
  coefmat <- cbind(
    `Estimate` = coef,
    `Std. Error` = asyse,
    `MCMC %` = est.pct,
    `z value` = zval,
    `Pr(>|z|)` = pval)

  rownames(coefmat) <- param_names(object)

  devtext <- "Deviance:"
  if (object$estimate!="MPLE" || !independence || object$reference != as.formula(~Bernoulli)) {
    if (pseudolikelihood) {
      devtext <- "Pseudo-deviance:"
      ans$message <- "\nWarning:  The standard errors are based on naive pseudolikelihood and are suspect.\n"
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
    dyads<- sum(as.rlebdm(object$constrained, object$constrained.obs, which="informative"))
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
