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














































































