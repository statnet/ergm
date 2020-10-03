#  File R/summary.ergm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################



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
#' @return The function [summary.ergm()] computes and returns a list of summary
#'   statistics of the fitted [ergm()] model given in `object`. Note that for
#'   backwards compatibility, it returns two coefficient tables: `$coefs` which
#'   does not contain the z-statistics and `$coefficeints` which does (and is
#'   therefore more similar to those returned by [stats::summary.lm()]).
#'
#'   The returned object is a list of class "ergm.summary" with the following
#'   elements:
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
#' \item{coefs, coefficients}{data frames with model parameters and associated statistics}
#' \item{asycov}{asymptotic covariance matrix}
#' \item{asyse}{asymptotic standard error matrix}
#' \item{offset, drop, estimate, iterations, mle.lik, null.lik}{
#' see documentation of the object returned by [ergm()]
#' }
#'
#' @seealso The model fitting function [ergm()], [print.ergm()], and
#'   [base::summary()]. Function [stats::coef()] will extract the data frame of
#'   coefficients with standard errors, t-statistics and p-values.
#'
#'
#'
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
  
  if(any(is.na(object$coef)) & !is.null(object$mplefit)){
     object$coef[is.na(object$coef)] <-
     object$mplefit$coef[is.na(object$coef)]
  }

  
  nodes<- network.size(object$network)

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
              control=object$control)
  
  ans$samplesize <- switch(object$estimate,
                           EGMME = NVL3(control$EGMME.main.method, switch(.,
                             `Gradient-Descent`=control$SA.phase3_n,
                             stop("Unknown estimation method. This is a bug."))),
                           MPLE = NA,
                           CD=,
                           MLE = NVL3(control$main.method, switch(.,
                             CD=control$MCMC.samplesize,
                             `Stochastic-Approximation`=,
                               MCMLE=control$MCMC.samplesize,
                             `Robbins-Monro`=control$RM.phase3n,
                             `Stepping`=control$Step.MCMC.samplesize,
                             stop("Unknown estimation method. This is a bug."))),
                           stop("Unknown estimate type. This is a bug.")
                           )
                              

  ans$iterations <- switch(object$estimate,
                           EGMME = NVL3(control$EGMME.main.method, switch(.,
                             `Gradient-Descent`=NA,
                             stop("Unknown estimation method. This is a bug."))),
                           MPLE = NA,
                           CD=control$CD.maxit,
                           MLE = NVL3(control$main.method, switch(.,
                               `Stochastic-Approximation`=NA,
                             MCMLE=paste(object$iterations, "out of", control$MCMLE.maxit),
                             CD=control$CD.maxit,
                             `Robbins-Monro`=NA,
                             `Stepping`=NA,
                             stop("Unknown estimation method. This is a bug."))),
                           stop("Unknown estimate type. This is a bug.")
                           )
  
  nodes<- network.size(object$network)
  dyads<- sum(as.rlebdm(object$constrained, object$constrained.obs, which="informative"))
  df <- length(object$coef)


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

  rdf <- dyads - df
  zval <- object$coef / asyse
  pval <- 2 * pnorm(q=abs(zval), lower.tail=FALSE)
  
  count <- 1
  coefmat <- cbind(
    `Estimate` = coef(object),
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
  null.lik<-try(logLikNull(object,...), silent=TRUE)

  ans$null.lik.0 <- is.na(null.lik)

  if(!inherits(mle.lik,"try-error")){

    ans$devtable <- matrix(c(if(is.na(null.lik)) 0 else -2*null.lik, -2*mle.lik,
                             c(dyads, rdf)), 2,2, dimnames=list(c("Null","Residual"),
                                                                c("Resid. Dev", "Resid. Df")))
    ans$devtext <- devtext
        
    ans$aic <- AIC(mle.lik)
    ans$bic <- BIC(mle.lik)
    ans$mle.lik <- ERRVL(mle.lik, NA)
    ans$null.lik <- ERRVL(null.lik, NA)
  }else ans$objname<-deparse(substitute(object))

  ans$coefs <- as.data.frame(coefmat)[,-3] # For backwards compatibility.
  ans$coefficients <- as.data.frame(coefmat)
  ans$asycov <- asycov
  ans$asyse <- asyse
  class(ans) <- "summary.ergm"
  ans
}

