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
    x <- matrix(x, ncol=length(x))
    x <- sweep(x, 2, center)
    cov <- robust.inverse(cov, ...)
    retval <- rowSums((x%*%cov) * x)
    names(retval) <- rownames(x)
    retval
}
