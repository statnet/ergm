#  File R/ergm.eta.R in package ergm, part of the Statnet suite of packages for
#  network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################

#' Operations to map curved [ergm()] parameters onto canonical parameters
#' 
#' This family of functions concerns mapping model parameter
#' \eqn{q}-vector \eqn{\theta} using to the canonical parameter
#' \eqn{p}-vector \eqn{\eta} based on the `etamap` object, usually
#' attached as the `$etamap` element of an [`ergm_model`] object.
#' 
#' These functions are mainly important in the case of curved exponential family
#' models, i.e., those in which the parameter of interest (\eqn{\theta}) is not a
#' linear function of the natural parameters (\eqn{\eta}) in the exponential-family
#' model. In non-curved models, we may assume without loss of generality that
#' \eqn{\eta(\theta)=\theta}.
#'
#' A succinct description of how eta(theta) is incorporated into an ERGM is
#' given by equation (5) of \insertCite{Hu07c;textual}{ergm}.  See \insertCite{HuHa06i;textual}{ergm} and
#' \insertCite{Hu07c;textual}{ergm} for further details about how eta and its derivatives are used
#' in the estimation process.
#' 
#' @param theta a curved model parameter \eqn{q}-vector
#' @param etamap the list of values that describes the \eqn{\theta\to\eta}{theta -> eta}
#'   mapping, usually attached as `$etamap` element of an [`ergm_model`]
#'   object. At this time, it is a list with the following elements:
#' \describe{
#' \item{`canonical`}{ a numeric vector whose `i`th entry specifies whether the `i`th component of theta is canonical (via non-negative integers) or curved (via zeroes)}
#' \item{`offsetmap`}{ a logical vector whose `i`th entry tells whether the ith coefficient of the canonical parameterization was "offset", i.e fixed}
#' \item{`offset`}{ a logical vector whose ith entry tells whether the ith model term was offset/fixed}
#' \item{`offsettheta`}{ a logical vector whose ith entry tells whether the ith curved theta coeffient was offset/fixed;}
#' \item{`curved`}{ a list with one component per curved EF term in
#' the model containing \describe{
#' \item{`from`}{ the indices of the curved theta parameter that are to be mapped from}
#' \item{`to`}{ the indices of the canonical eta parameters to be mapped to}
#' \item{`map`}{ the map provided by [`InitErgmTerm`]}
#' \item{`gradient`}{ the gradient function provided by [`InitErgmTerm`]}
#' \item{`cov`}{ optional additional covariates to be passed to the map and the gradient functions }
#' \item{`etalength`}{ the length of the eta vector}
#' }
#' }
#' }
#' @return For \code{ergm.eta}, the canonical \eqn{\eta} parameters as mapped
#'   from \eqn{\theta}.
#'
#' @seealso [`ergmTerm`]
#' @references \insertAllCited{}
#' @keywords internal
#' @name ergm.eta
NULL

#' @describeIn ergm.eta Evaluate \eqn{\eta(\theta):
#'   \mathbb{R}^q\to\mathbb{R}^p}{eta(theta): R^q -> R^p}.
#' @export ergm.eta
ergm.eta <- function(theta, etamap){
  .Call("ergm_eta_wrapper", as.numeric(theta), etamap, PACKAGE="ergm")
}

#' @describeIn ergm.eta Evaluate \eqn{\eta'(\theta):
#'   \mathbb{R}^q\to\mathbb{R}^{q\times p}}{deta(theta): R^q ->
#'   R^{q*p}}. If the gradient is only intended to be a multiplier for
#'   some vector or matrix, the more efficient \code{ergm.etagradmult}
#'   is recommended.
#' @export ergm.etagrad
ergm.etagrad <- function(theta, etamap){
  .Call("ergm_etagrad_wrapper", as.numeric(theta), etamap, PACKAGE="ergm")
}

#' @describeIn ergm.eta Evaluate \eqn{\eta'(\theta) V: \mathbb{R}^q \times \mathbb{R}^{p\times r}\to\mathbb{R}^{q\times
#'   r}}{deta(theta)\%*\%V: R^q * R^{p*r} -> R^{q*r}}.
#' @param v a \eqn{p}-vector the vector or a \eqn{p\times r}{p*r}
#'   matrix
#' @export ergm.etagradmult
ergm.etagradmult <- function(theta, v, etamap){
  storage.mode(v) <- "double"
  .Call("ergm_etagradmult_wrapper", as.numeric(theta), v,  etamap, PACKAGE="ergm")
}

#' @describeIn ergm.eta A convenience function to evaluate
#'   \eqn{\{\eta'(\theta) S^{\top}\}^{\top} = S
#'   \eta'(\theta)^{\top}}{t(deta(theta)\%*\%t(S)) =
#'   S\%*\%t(deta(theta))}.
#' @param s a \eqn{r\times p} matrix
#' @export ergm.etagradmultt
ergm.etagradmultt <- function(theta, s, etamap) {
  storage.mode(s) <- "double"
  ## TODO: Write a C-side version to avoid transposes.
  t(.Call("ergm_etagradmult_wrapper", as.numeric(theta), t(s),  etamap, PACKAGE = "ergm"))
}
