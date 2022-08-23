#  File R/ergm.eta.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
##############################################################################
# The <ergm.eta> function calculates and returns eta, mapped from
# theta using the etamap object created by <ergm.etamap>.
#
# --PARAMETERS--
#   theta :  the curved model parameters  
#   etamap:  the list of values that constitutes the theta-> eta mapping
#            and is returned by <ergm.etamap>
#
# --RETURNED--
#   eta:  the canonical eta parameters as mapped from theta
#
###############################################################################



#' Operations to map curved [ergm()] parameters onto canonical parameters
#' 
#' The \code{ergm.eta} function calculates and returns eta, mapped
#' from theta using the `etamap` object, usually attached as the
#' `$etamap` element of an [`ergm_model`] object.
#' 
#' These functions are mainly important in the case of curved exponential family
#' models, i.e., those in which the parameter of interest (theta) is not a
#' linear function of the natural parameters (eta) in the exponential-family
#' model. In non-curved models, we may assume without loss of generality that
#' eta(theta)=theta.
#' 
#' A succinct description of how eta(theta) is incorporated into an ERGM is
#' given by equation (5) of Hunter (2007).  See Hunter and Handcock (2006) and
#' Hunter (2007) for further details about how eta and its derivatives are used
#' in the estimation process.
#' 
#' @param theta the curved model parameters
#' @param etamap the list of values that describes the theta -> eta
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
#' @return For \code{ergm.eta}, the canonical eta parameters as mapped
#'   from theta.
#'
#' @seealso \code{\link{ergmTerm}}
#' @references \itemize{ \item Hunter, D. R. and M. S. Handcock
#'   (2006).  Inference in curved exponential family models for
#'   networks. \emph{Journal of Computational and Graphical
#'   Statistics}, 15: 565--583.
#' 
#' \item Hunter, D. R. (2007). Curved exponential family models for social
#' networks. \emph{Social Networks}, 29: 216--230.  }
#' @keywords internal
#' @export ergm.eta
ergm.eta <- function(theta, etamap){
  .Call("ergm_eta_wrapper", as.numeric(theta), etamap, PACKAGE="ergm")
}
