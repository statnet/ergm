#  File R/is.inCH.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2015 Statnet Commons
#######################################################################
###############################################################################
# The <is.inCH> function determines whether a vector p is in the convex hull
# of the vectors M
#
# --PARAMETERS--
#   p:  a vector of length n
#   M:  a q by n matrix 
#
# --RETURNED--
#   x: TRUE if p is in the CH of the points M
#      FALSE othewise
#
###############################################################################

########
# The n-vector p is in the convex hull of the n-vectors in the
# Rxn matrix M iff the maximum value of z_0 + z'p equals zero for all vectors z
# satisfying z_0 + z'M \le 0.  Thus, the value returned by is.inCH depends on the
# solution of a linear program, for which the lpSolve package and its function
# lp is needed.  Letting q=(1 p')' and L = (1 M), the program is:
#
#  maximize z'q  
#  subject to z'L <= 0 
#
#  Notice that, if there exists a z that makes z'q positive, then there is
#  no maximizer since in that case Kz gives a larger value whenever K>1.  
#  For this reason, we add one additional constraint, namely, z'q <= 1.
#
#  To put this all in "standard form", we let z=a-b, where a, b, are nonnegative.
#  If we then write x=(a' b')', we obtain a new linear program:  
#
#  Minimize x'(-q' q')'
#  subject to x'(q' -q')' <= 1 and x'(L -L) <= 0 and x >= 0
#  ...and if the minimum is strictly negative, return FALSE because the point
#  is not in the CH in that case.

## Note: p can be a matrix. In that case, every row of p is checked.

is.inCH <- function(p, M, verbose=FALSE, ...) { # Pass extra arguments directly to LP solver

  if(is.null(dim(p))) p <- rbind(p)

  if (!is.matrix(M)) 
    stop("Second argument must be a matrix.")
  if (ncol(p) != ncol(M)) 
    stop("Number of columns in matrix (2nd argument) is not equal to dimension ",
         "of first argument.")

  if(nrow(M)==1){
    for(i in seq_len(nrow(p))){
      if(!isTRUE(all.equal(p[i,], M, check.attributes = FALSE))) return(FALSE)
    }
    return(TRUE)
  }

  ##
  ## NOTE: PCA code has been moved to .Hummel.steplength().
  ##

  L = cbind(1, M)

  for(i in seq_len(nrow(p))){
   q = c(1, p[i,]) 
############################################
# USE lp FUNCTION FROM lpSolve PACKAGE:
   ans <- lp(objective.in = c(-q, q),
             const.mat = rbind( c(q, -q), cbind(L, -L)),
             const.dir = "<=",
             const.rhs = c(1, rep(0, NROW(L))), 
             ...
             )
   if(ans$objval!=0){
    if(verbose) cat(sprintf("is.inCH: iter= %d, outside hull.\n",i))
    return(FALSE)  #if the min is not zero, the point p[i,] is not in the CH of the points M
   }
  }
  if(verbose) cat(sprintf("is.inCH: iter= %d, inside hull.\n",i))
  return(TRUE) # If all points passed the test, return TRUE.

##############################################
## USE solveLP FUNCTION FROM linprog PACKAGE (deprecated)
## From help for function 'solveLP' in package 'linprog':
##     Minimizes (or maximizes) c'x, subject to A x <= b and x >= 0.
#  ans <- solveLP (cvec = c(-q, q),
#                  bvec = c(1, rep(0, NROW(L))),
#                  Amat= rbind( c(q, -q), cbind(L, -L)),
#                  ...
#                  )
#  if(ans$opt==0)return(TRUE)  #if the min is zero, the point p is in the CH of the points M
#  else return(FALSE)              

### OLD CODE USING Rglpk PACKAGE (deprecated)
#  R = NROW(M)
#  C=length(p)+1
#	ans <- Rglpk_solve_LP(obj=c(p,-1), 
#	                      mat=cbind(rbind(p,M),-1),
#	                      dir=as.vector(rep("<=",R+1)), 
#	                      rhs=as.vector(c(1,rep(0,R))),
#	                      max=TRUE, 
#	                      bounds=list(lower=list(ind=1:C,val=rep(-Inf,C))))
#  if(ans$optimum==0)return(TRUE)  #if the max is zero, the point p is in the CH of the points M
#  else return(FALSE)

}
