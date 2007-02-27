ergm.mahalanobis <- function(x, center, cov, inverted=FALSE, ...)
{
    x <- matrix(x, ncol=length(x))
    x <- sweep(x, 2, center)
    cov <- robust.inverse(cov, ...)
    retval <- rowSums((x%*%cov) * x)
    names(retval) <- rownames(x)
    retval
}
