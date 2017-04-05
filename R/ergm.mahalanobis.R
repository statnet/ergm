#  File ergm/R/ergm.mahalanobis.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
################################################################
# The <ergm.mahalanobis> function computes and returns the
# mahalanobis distance
################################################################

ergm.mahalanobis <- function(x, center, cov, inverted=FALSE, ...)
{
    x <- matrix(x, ncol=length(x))
    x <- sweep(x, 2, center)
    cov <- robust.inverse(cov, ...)
    retval <- rowSums((x%*%cov) * x)
    names(retval) <- rownames(x)
    retval
}
