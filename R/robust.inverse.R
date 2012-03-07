#  File ergm/R/robust.inverse.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
###############################################################################
# The <robust.inverse> function attempts to return the inverse of a matrix H
# by either direct means (via R's <solve>) or via computation from the
# singular decomposition of H
###############################################################################

"robust.inverse" <- function (H, tol = sqrt(.Machine$double.eps)) 
{
    iH <- try(solve(H), silent=TRUE)
    if(inherits(iH,"try-error")){
     if (length(dim(H)) > 2 || !(is.numeric(H) || is.complex(H))) 
        stop("H must be a numeric or complex matrix")
     if (!is.matrix(H)) 
        H <- as.matrix(H)
     Hsvd <- try(svd(H), silent=TRUE)
     if(inherits(Hsvd,"try-error")){
        warning("ergm did not compute all the standard errrors.")
        return(H)
     }
     if (is.complex(H)) 
        Hsvd$u <- Conj(Hsvd$u)
     Positive <- Hsvd$d > max(tol * Hsvd$d[1], 0)
     if (all(Positive)) 
        Hsvd$v %*% (1/Hsvd$d * t(Hsvd$u))
     else if (!any(Positive)) 
        array(0, dim(H)[2:1])
     else Hsvd$v[, Positive, drop = FALSE] %*% ((1/Hsvd$d[Positive]) * 
         t(Hsvd$u[, Positive, drop = FALSE]))
    }else{
     iH
    }
}
