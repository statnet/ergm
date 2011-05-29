##############################################################################
# The <ergm.theta.sample> function calculates and returns eta, mapped from
# theta using the etamap object created by <ergm.etamap>.
#
# --PARAMETERS--
#   theta :  the curved model parameters  
#   etamap:  the list of values that constitutes the theta-> eta mapping
#            and is returned by <ergm.etamap>
#   sample:  MCMC sample statistics returned by <ergm.estimate>
#
# --RETURNED--
#   sample:  MCMC sample statistics from curved model (conditional on MLE)
#
###############################################################################

ergm.theta.sample <- function(theta, etamap, sample) {
  s <- sample
  if(length(etamap$curved)>0) {
    eta <- rep(0,etamap$etalength)
    ec <- etamap$canonical
    eta[ec[ec>0]] <- theta[ec>0]
    delete <- rep(FALSE,etamap$etalength)
    delete.names <- rep(FALSE,length(theta))
    for(i in 1:length(etamap$curved)) {
      cm <- etamap$curved[[i]]
      loc <- min(cm$to)
      delete[cm$to[-which.min(cm$to)]] <- TRUE
      eta[cm$to] <- cm$map(theta[cm$from],length(cm$to),cm$cov)
      s[,loc] <- sample[,cm$to] %*% eta[cm$to]
#  MSH note: The next line only works for canonical parameter is first
#  It is ok for the current curved stats
      s[,loc] <- s[,loc] / theta[cm$from][1]
      loc <- min(cm$from)
      delete.names[cm$from[-which.min(cm$from)]] <- TRUE
    }
    s <- s[,!delete]
    colnames(s) <- names(theta)[!delete.names]
  }
  s
}
