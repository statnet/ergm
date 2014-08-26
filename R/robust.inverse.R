#  File R/robust.inverse.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
###############################################################################
# The <robust.inverse> function attempts to return the inverse of a matrix H
# by either direct means (via R's <solve>) or via computation from the
# singular decomposition of H
#
# --PARAMETERS--
#   H  : a matrix, presumably a Hessian
#   tol: the tolerance, used for determining which of the singular values
#        from the decompostion are postive; default=sqrt(.Machine$double.eps)
#
# --RETURNED--
#   H : the original matrix, if H could not be inverted directly and the
#       singular decomposition of H incurred an error
#   0 : a p by q matrix of 0's, where H is q by p, if none of the singular
#       values from the decomposition are positive
#   iH: the inverse of H, computed directly if possible, otherwise computed
#       from the components of the singular decomposition
#
###############################################################################

.robust.inverse <- function (H, tol = sqrt(.Machine$double.eps)) 
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
