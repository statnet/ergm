#  File ergm/R/ergm.eta.R
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
ergm.eta <- function(theta, etamap) {
# This function maps theta to eta based on the etamap object created
# by ergm.etamap.
#eta.cov added by CTB on 1/28/06
  eta <- rep(0,etamap$etalength)
  ec <- etamap$canonical
  eta[ec[ec>0]] <- theta[ec>0]
  if(length(etamap$curved)>0) {
    for(i in 1:length(etamap$curved)) {
      cm <- etamap$curved[[i]]
      eta[cm$to] <- cm$map(theta[cm$from],length(cm$to),cm$cov)
    }
  }
  eta[is.infinite(eta)] <- -10000
  eta
}














































































