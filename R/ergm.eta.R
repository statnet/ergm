#  File ergm/R/ergm.eta.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
##############################################################################
# The <ergm.eta> function calculates and returns eta, mapped from
# theta using the etamap object created by <ergm.etamap>.
###############################################################################

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
  eta[is.infinite(eta)] <- sign(eta[is.infinite(eta)])*1000
  eta
}














































































