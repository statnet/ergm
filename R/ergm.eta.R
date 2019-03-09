#  File R/ergm.eta.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################
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
#' The \code{ergm.eta} function calculates and returns eta, mapped from theta
#' using the etamap object created by \code{ergm.etamap}.
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
#' @param etamap the list of values that constitutes the theta-> eta
#'   mapping and is returned by \code{ergm.etamap}
#' @return For \code{ergm.eta}, the canonical eta parameters as mapped
#'   from theta.
#' @seealso \code{\link{ergm-terms}}
#' @references \itemize{ \item Hunter, D. R. and M. S. Handcock
#'   (2006).  Inference in curved exponential family models for
#'   networks. \emph{Journal of Computational and Graphical
#'   Statistics}, 15: 565--583.
#' 
#' \item Hunter, D. R. (2007). Curved exponential family models for social
#' networks. \emph{Social Networks}, 29: 216--230.  }
#' @keywords internal
#' @export ergm.eta
ergm.eta <- function(theta, etamap) {
  eta <- rep(0,etamap$etalength)
  ec <- etamap$canonical
  eta[ec[ec>0]] <- theta[ec>0]
  if(length(etamap$curved)>0) {
    for(i in 1:length(etamap$curved)) {
      cm <- etamap$curved[[i]]
      eta[cm$to] <- cm$map(theta[cm$from],length(cm$to),cm$cov)
    }
  }
  eta
}














































































