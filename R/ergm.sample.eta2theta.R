#  File R/ergm.sample.eta2theta.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
##############################################################################
# The <ergm.theta.sample> function calculates and returns eta, mapped from
# theta using the etamap object created by <ergm.etamap>.
#
# --PARAMETERS--
#   sample:  MCMC sample statistics returned by <ergm>
#   coef :  the curved model parameters  
#   etamap:  the list of values that constitutes the theta-> eta mapping
#            and is returned by <ergm.etamap>
#
# --RETURNED--
#   sample:  MCMC sample statistics from curved model (conditional on MLE)
#
###############################################################################


ergm.sample.eta2theta <- function(sample, coef, etamap) {

  which.main <- list(gwdsp=1,gwesp=1,gwnsp=1,
                     altkstar=1,
                     gwb1degree=1,gwb2degree=1,
                     gwdegree=1,gwidegree=1,gwodegree=1)

  
  s <- sample
  if(length(etamap$curved)>0) {
    eta <- rep(0,etamap$etalength)
    ec <- etamap$canonical
    eta[ec[ec>0]] <- coef[ec>0]
    delete <- rep(FALSE,etamap$etalength)
    delete.names <- rep(FALSE,length(coef))
    for(cm in etamap$curved) {
      # Figure out which, if any, element of coef is the "main" one in the curved expression.
      main.name<-names(coef)[cm$from][names(coef)[cm$from] %in% names(which.main)]
      if(!length(main.name)){
        warning("Curved term with parameters ",paste(names(coef)[cm$from],collapse=" and "), " is not in the whitelist and was not mapped.")
        break
      }
      main.pos<-which.main[[main.name]]

      loc <- min(cm$to)
      delete[cm$to[-which.min(cm$to)]] <- TRUE
      eta[cm$to] <- cm$map(coef[cm$from],length(cm$to),cm$cov)
      s[,loc] <- sample[,cm$to] %*% eta[cm$to]
#  MSH note: The next line only works for canonical parameter is first
#  It is ok for the current curved stats
      s[,loc] <- s[,loc] / coef[cm$from][main.pos]
      delete.names[cm$from[-main.pos]] <- TRUE
    }
    s <- s[,!delete]
    colnames(s) <- names(coef)[!delete.names]
  }
  s
}
