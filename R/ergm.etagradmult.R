#  File R/ergm.etagradmult.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2015 Statnet Commons
#######################################################################
##############################################################################
# The <ergm.etagradmult> function calculates and returns the product of the
# gradient of eta with a vector v
#
# --PARAMETERS--
#   theta :  the curved model parameters
#   v     :  a vector of the same length as the vector of mapped eta parameters
#   etamap:  the list constituting the theta-> eta mapping returned by <ergm.etamap>
#
# --RETURNED--
#   ans: the vector that is the product of the gradient of eta and v; infinite
#        values are replaced by (+-)10000
#
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
  ans
}













































































