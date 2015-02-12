#  File R/is.inCH.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
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

is.inCH <- function(p, M) {
  p <- as.vector(p)
  if (!is.matrix(M)) 
    stop("Second argument must be a matrix.")
  if (length(p) != NCOL(M)) 
    stop("Number of columns in matrix (2nd argument) is not equal to dimension ",
         "of first argument.")
  R = NROW(M)
  C=length(p)+1
  
	ans <- Rglpk_solve_LP(obj=c(p,-1), mat=cbind(rbind(p,M),-1),
                              dir=as.vector(rep("<=",R+1)), rhs=as.vector(c(1,rep(0,R))),
                              max=TRUE, bounds=list(lower=list(ind=1:C,val=rep(-Inf,C))))
  
  if(ans$optimum==0){x<-TRUE}  #if the max is zero, the point q is in the CH of the points p
  else{x<-FALSE}              #if the max is strictly positive, q is not in the CH of p
  x
}




