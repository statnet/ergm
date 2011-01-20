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
#   eta:  the canonical eta parameters as mapped from theta; infinite parameters
#         are replaced by -10000
#
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
  eta[is.infinite(eta)] <- -10000
  eta
}














































































