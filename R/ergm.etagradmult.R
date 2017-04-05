#  File ergm/R/ergm.etagradmult.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
##############################################################################
# The <ergm.etagradmult> function calculates and returns the product of the
# gradient of eta with a vector v
#################################################################################

ergm.etagradmult <- function(theta, v, etamap) {
  v <- as.matrix(v)
  ans <- matrix(0, length(theta), dim(v)[2])
  if(dim(v)[1] != etamap$etalength)
    stop("Non-conforming matrix multiply: grad(eta) %*% v.\n",
         "grad(eta) has ", etamap$etalength, " columns ",
         "and v has ", dim(v)[1], " rows.")
  ec <- etamap$canonical
# Set gradient for canonical parameters to the identity matrix
  ans[ec>0,] <- v[ec[ec>0],]
  if(length(etamap$curved)>0) {
    for(i in 1:length(etamap$curved)) {
      cm <- etamap$curved[[i]]
      #cov added by CTB on 1/28/06
      ans[cm$from,] <- cm$gradient(theta[cm$from], length(cm$to), cm$cov)%*%v[cm$to,]  
    }
  }
  ans[is.infinite(ans)] <- 10000*sign(ans)[is.infinite(ans)]
  ans
}













































































