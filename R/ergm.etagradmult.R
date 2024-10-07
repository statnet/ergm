#  File R/ergm.etagradmult.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################
##############################################################################
# The <ergm.etagradmult> function calculates and returns the product of the
# gradient of eta with a vector v
#
# --PARAMETERS--
#   theta :  the curved model parameters
#   v     :  a vector of the same length as the vector of mapped eta parameters
#   etamap:  the list constituting the theta-> eta mapping returned by <ergm.etamap>
#
# --RETURNED--
#   ans: the vector that is the product of the gradient of eta and v; infinite
#        values are replaced by (+-)10000
#
#################################################################################

#' @rdname ergm.eta
#' @description The \code{ergm.etagradmult} function calculates and
#'   returns the product of the gradient of eta with a vector `v`.
#' @param v a vector of the same length as the vector of mapped eta
#'   parameters
#' @return For \code{ergm.etagradmult}, the vector that is the product
#'   of the gradient of eta and \code{v}.
#' @export ergm.etagradmult
ergm.etagradmult <- function(theta, v, etamap){
  storage.mode(v) <- "double"
  .Call("ergm_etagradmult_wrapper", as.numeric(theta), v,  etamap, PACKAGE="ergm")
}
