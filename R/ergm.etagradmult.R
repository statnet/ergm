ergm.etagradmult <- function(theta, v, etamap) {
# This function returns g %*% v, where g is the gradient of eta(theta)
# based on the etamap object created by ergm.etamap.
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
      ans[cm$from,] <- cm$gradient(theta[cm$from], length(cm$to), cm$cov)%*%v[cm$to,]  #cov added by CTB on 1/28/06
    }
  }
  ans[is.infinite(ans)] <- 10000*sign(ans)[is.infinite(ans)]
  ans
}













































































