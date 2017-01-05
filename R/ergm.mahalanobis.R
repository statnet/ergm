#  File R/ergm.mahalanobis.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
################################################################
# The <ergm.mahalanobis> function computes and returns the
# mahalanobis distance
#
# --PARAMETERS--
#   x     : a random vector of length n
#   center: the mean values of 'x'
#   cov   : the covariance matrix of 'x'  
#
# --IGNORED--
#   inverted: ?? 
#
# --RETURNED--
#   retval: the mahalanobis distance
#
################################################################

ergm.mahalanobis <- function(x, center, cov, inverted=FALSE, ...)
{
  .Deprecated(msg = 'ergm.mahalanobis will be removed in future versions of ergm. See stats::mahalanobis as an alternative')
    x <- matrix(x, ncol=length(x))
    x <- sweep(x, 2, center)
    cov <- robust.inverse(cov, ...)
    retval <- rowSums((x%*%cov) * x)
    names(retval) <- rownames(x)
    retval
}
