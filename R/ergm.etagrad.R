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















































































