#  File R/ergm.etagrad.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################
################################################################################
# The <ergm.etagrad> function caculates and returns the gradient of eta
# mapped from theta using the etamap object created by <ergm.etamap>. If the
# gradient is only intended to be a multiplier for some vector, the more
# efficient <ergm.etagradmult> is recommended.
#
# --PARAMETERS--
#   theta :  the vector of curved model parameters 
#   etamap:  the list constituting the theta-> eta mapping that is returned by
#            <ergm.etamap>
#
# --RETURNED--
#   etagrad: a matrix of the gradient of eta 
#
################################################################################

#' @rdname ergm.eta
#' @description The \code{ergm.etagrad} function caculates and returns
#'   the gradient of eta mapped from theta using the etamap object
#'   created by \code{ergm.etamap}. If the gradient is only intended
#'   to be a multiplier for some vector, the more efficient
#'   \code{ergm.etagradmult} is recommended.
#' @return For \code{ergm.etagrad}, a matrix of the gradient of eta
#'   with respect to theta.
#' @export ergm.etagrad
ergm.etagrad <- function(theta, etamap) {
  etagrad <- matrix(0, length(theta), etamap$etalength)
  ec <- etamap$canonical
# Set gradient for canonical parameters to the identity matrix
  etagrad[ec>0, ec[ec>0]] <- diag(sum(ec>0))
  if(length(etamap$curved)>0) {
    for(i in 1:length(etamap$curved)) {
      cm <- etamap$curved[[i]]
      etagrad[cm$from,cm$to] <- cm$gradient(theta[cm$from], length(cm$to), cm$cov)  #Added by CTB on 1/28/06
    }
  }
  etagrad
}















































































