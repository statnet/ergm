#  File ergm/R/ergm.etagrad.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
################################################################################
# The <ergm.etagrad> function caculates and returns the gradient of eta
# mapped from theta using the etamap object created by <ergm.etamap>. If the
# gradient is only intended to be a multiplier for some vector, the more
# efficient <ergm.etagradmult> is recommended.
################################################################################

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
  etagrad[is.infinite(etagrad)] <- 10000*sign(etagrad)[is.infinite(etagrad)]
  etagrad
}















































































