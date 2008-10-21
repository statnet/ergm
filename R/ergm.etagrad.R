#  File ergm/R/ergm.etagrad.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of Washington
#                David R. Hunter, Penn State University
#                Carter T. Butts, University of California - Irvine
#                Steven M. Goodreau, University of Washington
#                Martina Morris, University of Washington
# Copyright 2007 The statnet Development Team
######################################################################
ergm.etagrad <- function(theta, etamap) {
# This function maps theta to the gradient of eta(theta) based on the etamap
#  object created by ergm.etamap.
# Note that if you only need the gradient in order to multiply it by a
# vector, it is much more efficient to use the ergm.etagradmult function.
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















































































