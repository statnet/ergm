ergm.checkdegeneracy <- function(statsmatrix, statsmatrix.miss=NULL, verbose=FALSE) {
#
 degen <- FALSE
 novar <- apply(statsmatrix,2,var)<1e-6
 if(all(novar)){
  if(verbose){
    warning("All the MCMC sample statistics are the same.\n")
    print(apply(statsmatrix,2,summary.statsmatrix.ergm),scipen=6)
  }
  degen <- TRUE
 }
# 
#  ########### from ergm.estimate:
#   # First, perform simple check to see if observed statistics
#   # (all zeros) lie outside the range of the MCMC sample statistics.
#  outofbounds <- (apply(statsmatrix,2,max)*
#                  apply(statsmatrix,2,min)) > epsilon
#  if (any(outofbounds)){
#    ergm.marquardt()
#  }
#
#  ############# from ergm.statseval:  CHECK FOR NO VARIANCE
#  # below is what was done; not necessarily something to keep
#  novar=rep(FALSE,length(theta0)) # This is a hack
#  theta0[!novar] <- l$coef
#  l$coef <- theta0
#  theta0[!novar] <- l$MCMCtheta
#  covar <- 0*diag(length(l$coef)) # initialize to zero matrix
#  covar[!novar,!novar] <- l$covar
#
#
#  ############## from ergm.mainfitloop:
#  if(sum(z$newg) > 12000){
#    cat("Returned network is too full. Retaining the original.\n")
#    newnetwork <- g
#  }
degen
}
