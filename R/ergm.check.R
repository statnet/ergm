#  File R/ergm.check.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################
#============================================================================
# This file contains the following 4 functions for checking ergm arguments
#          <ergm.checkargs>           <ergm.checkdirected>
#          <ergm.checkbipartite>      <ergm.checkdegeneracy>
#============================================================================





######################################################################################
# The <ergm.checkargs> function ensures for the <InitErgm.X> function that the
# term X:
#   1) has an appropiate number of arguments
#   2) has correct argument types if arguments where provieded
#   3) has default values assigned for non-required arguments
# by halting execution if either of the first 2 criteria are not met
#
# --PARAMETERS--
#  fname        : the name of the model term as a character string
#  arglist      : the list of arguments for term X
#  varnames     : the vector of names of the possible arguments for term X;
#                 default=NULL 
#  vartypes     : the vector of types of the possible arguments for term X;
#                 default=NULL 
#  defaultvalues: the list of default values for the possible arguments of term X;
#                 default=list()
#  required     : the logical vector of whether each possible argument is required;
#                 default=NULL
#
# --RETURNED--
#   a vector containing the following 2 variables:
#     .conflicts.OK: always TRUE
#     out          : a list of the values for each possible argument of term X;
#                    user provided values are used when given, default values otherwise.
#
######################################################################################

#' @include ergm-deprecated.R

#################################################################################
# The <ergm.checkdegeneracy> function checks for degeneracy by looking for
# variability in the stats matrix 
#
# --PARAMETERS--
#   statsmatrix :  the matrix of summary sample statistics
#   verbose     :  whether the degeneracy warning should be printed (T or F)
#
#
# --IGNORED PARAMETERS--
#   statsmatrix.obs : the matrix of summary sample statistics from the
#                      observation process network; default=NULL
#
#
# --RETURNED--
#   degen: whether the ergm model is degenerate, in the sense that
#          there was no varibility in the statsmatrix (T or F)
#
#################################################################################

ergm.checkdegeneracy <- function(statsmatrix, statsmatrix.obs=NULL, verbose=FALSE) {
 degen <- FALSE
 novar <- apply(statsmatrix,2,stats::var)<.Machine$double.eps^0.5
 if(all(novar)){
  if(verbose){
    warning("All the MCMC sample statistics are the same.\n", call.=FALSE)
    print(apply(statsmatrix,2,summary.statsmatrix.ergm),scipen=6)
  }
  degen <- TRUE
 }
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
#  novar=rep(FALSE,length(init)) # This is a hack
#  init[!novar] <- l$coef
#  l$coef <- init
#  init[!novar] <- l$MCMCtheta
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
